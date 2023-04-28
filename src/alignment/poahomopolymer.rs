
use petgraph::{graph::{NodeIndex}, visit::Topo, Directed, Graph, Incoming, Outgoing};

use crate::misc::HomopolymerCell;

pub const MIN_ISIZE: isize = -858_993_459;
pub const MAX_USIZE: usize = 858_993_459;

pub struct Poa {
    pub poa_graph: Graph<HomopolymerCell, i32, Directed, usize>,
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
}

#[derive(Clone)]
pub struct SimpleMatrixCell {
    score: isize,
    prev: usize,
    freq: usize,
    back: char,
}

impl Poa {
    pub fn initialize (start_query: &Vec<HomopolymerCell>, match_score: i32, mismatch_score: i32, gap_open_score: i32) -> Self {
        // allocate space for the graph
        let mut graph: Graph<HomopolymerCell, i32, Directed, usize> = Graph::with_capacity(start_query.len(), start_query.len() - 1);
        // add nodes and edges to the graph
        let mut prev: NodeIndex<usize> = graph.add_node(start_query[0].clone());
        let mut node: NodeIndex<usize>;
        for base in start_query.iter().skip(1) {
            node = graph.add_node(base.clone());
            graph.add_edge(prev, node, 1);
            prev = node;
        }
        Poa {poa_graph: graph, match_score: match_score, mismatch_score: mismatch_score, gap_open_score: gap_open_score}
    }
    pub fn add_to_poa (&mut self, query: &Vec<HomopolymerCell>) {
        // get topological ordering of the graph
        let mut topo = Topo::new(&self.poa_graph);
        let mut topo_indices = Vec::new();
        while let Some(node) = topo.next(&self.poa_graph) {
            topo_indices.push(node);
        }
        let align_vec = self.get_alignment_simple(query);
        let mut graph_index = topo_indices.len();
        let mut query_index = query.len();
        let mut current_node: Option<NodeIndex<usize>> = None;
        let mut prev_node: Option<NodeIndex<usize>>;
        let mut prev_node_require_update = false;
        for alignment in align_vec {
            // find the current processing node in topoindices
            match alignment.0 as char {
                'i' => {
                    graph_index = alignment.1;
                },
                'm' => {
                    prev_node = current_node;
                    current_node = Some(topo_indices[graph_index - 1]);
                    let _temp_base = self.poa_graph.node_weight(current_node.unwrap()).unwrap().base;
                    // modify this
                    //*self.poa_graph.node_weight_mut(current_node.unwrap()).unwrap() = HomopolymerCell::new(_temp_base, alignment.2);
                    if prev_node_require_update {
                        match self.poa_graph.find_edge(current_node.unwrap(), prev_node.unwrap()) {
                            Some(edge) => {
                                *self.poa_graph.edge_weight_mut(edge).unwrap() += 1;
                            }
                            None => {
                                self.poa_graph.add_edge(current_node.unwrap(), prev_node.unwrap(), 1);
                            }
                        }
                    }
                    query_index -= 1;
                    graph_index = alignment.1;
                    prev_node_require_update = true;
                },
                's' => {
                    prev_node = current_node;
                    current_node = Some(self.poa_graph.add_node(query[query_index - 1].clone()));
                    match prev_node {
                        Some(x) => {
                            self.poa_graph.add_edge(current_node.unwrap(), x, 1);
                        },
                        None => {},
                    }
                    prev_node_require_update = true;
                    query_index -= 1;
                    graph_index = alignment.1;
                }
                'd' => {
                    // make a new node
                    prev_node = current_node;
                    current_node = Some(self.poa_graph.add_node(query[query_index - 1].clone()));
                    match prev_node {
                        Some(x) => {
                            // connect the new node to the current one
                            self.poa_graph.add_edge(current_node.unwrap(), x, 1);
                        },
                        None => {},
                    }
                    prev_node_require_update = true;
                    query_index -= 1;
                },
                _ => (),
            }
        }
    }

    pub fn get_alignment_simple (&self, query: &Vec<HomopolymerCell>) -> Vec<(u8, usize, usize)> {
        // get topological ordering of the graph
        let mut topo = Topo::new(&self.poa_graph);
        let mut topo_indices = Vec::new();
        while let Some(node) = topo.next(&self.poa_graph) {
            topo_indices.push(node);
        }
        // initialize the weight matrix with zeros
        let mut poa_matrix: Vec<Vec<SimpleMatrixCell>> = vec![vec![SimpleMatrixCell { score: (0), prev: (0), freq: (0), back: ('0') }; query.len() + 1]; topo_indices.len() + 1];
        // fill the first row and column
        let mut graph_index = 1;
        let mut query_index = 1;
        // the first column
        for graph_node in &topo_indices {
            if graph_index == 1 {
                let temp_freq = self.poa_graph.raw_nodes()[graph_node.index()].weight.frequency;
                let temp_value = (self.gap_open_score * temp_freq as i32) as isize;
                poa_matrix[graph_index][0] = SimpleMatrixCell { score: (temp_value), prev: (0), freq: (temp_freq), back: ('i') };
            }
            else if graph_index > 1 {
                let mut max_score = MIN_ISIZE;
                let mut max_position = 0;
                let temp_freq = self.poa_graph.raw_nodes()[graph_node.index()].weight.frequency;
                // go through all the incoming nodes to this one and find the one with the max score
                let prevs: Vec<NodeIndex<usize>> = self.poa_graph.neighbors_directed(*graph_node, Incoming).collect();
                for prev_node in &prevs {
                    // find the index of the node in topo list
                    let position = topo_indices.iter().position(|r| r == prev_node).unwrap();
                    // get the score and add gap_extend
                    let temp_score = poa_matrix[position + 1][0].score + (self.gap_open_score * temp_freq as i32) as isize;
                    //save the max score and position
                    if temp_score > max_score {
                        max_score = temp_score;
                        max_position = position + 1;
                    }
                }
                poa_matrix[graph_index][0] = SimpleMatrixCell { score: (max_score), prev: (max_position), freq: (temp_freq), back: ('i') };
            }
            graph_index += 1;
        }
        // the first row
        for base in query {
            if query_index == 1 {
                let temp_value = (self.gap_open_score * base.frequency as i32) as isize;
                poa_matrix[0][query_index] = SimpleMatrixCell { score: (temp_value), prev: (0), freq: (base.frequency), back: ('d') };
            }
            else if query_index > 1 {
                let temp_value = poa_matrix[0][query_index - 1].score + (self.gap_open_score * base.frequency as i32) as isize;
                poa_matrix[0][query_index] = SimpleMatrixCell { score: (temp_value), prev: (0), freq: (base.frequency), back: ('d') };
            }
            query_index += 1;
        }
        // fill the rest
        for graph_index in 1..topo_indices.len() + 1 {
            for query_index in 1..query.len() + 1 {
                let temp_freq_query = query[query_index - 1].frequency;
                let temp_freq_graph = self.poa_graph.raw_nodes()[topo_indices[graph_index - 1].index()].weight.frequency;
                let temp_freq_match = (temp_freq_graph + temp_freq_graph) / 2 + 1;

                let temp_del_score = poa_matrix[graph_index][query_index - 1].score + (self.gap_open_score * temp_freq_query as i32) as isize;
                let mut temp_ins_score = MIN_ISIZE;
                let mut temp_ins_position = 0;
                let mut temp_match_score = MIN_ISIZE;
                let mut temp_match_position = 0;
                let prevs: Vec<NodeIndex<usize>> = self.poa_graph.neighbors_directed(topo_indices[graph_index - 1], Incoming).collect();
                for prev_node in &prevs {
                    // find the index of the node in topo list
                    let position = topo_indices.iter().position(|r| r == prev_node).unwrap();
                    //println!("position {} == prevnodeindex {}", position, prev_node.index());
                    // highest insert score and location
                    // get the score and add gap_extend
                    let temp_score = poa_matrix[position + 1][query_index].score + (self.gap_open_score * temp_freq_graph as i32) as isize;
                    // save the max score and position
                    if temp_score > temp_ins_score {
                        temp_ins_score = temp_score;
                        temp_ins_position = position + 1;
                    }
                    // highest match score and location
                    let temp_score = poa_matrix[position + 1][query_index - 1].score + (self.match_score * temp_freq_match as i32) as isize;
                    // save the max score and position
                    if temp_score > temp_match_score {
                        //println!("match score : {} ", temp_score - self.match_score as isize);
                        temp_match_score = temp_score;
                        temp_match_position = position + 1;
                    }
                }
                //println!("match score : {} ", temp_match_score - self.match_score as isize);
                //println!("test pos = {} graph index - 1 == {}", test_position, graph_index - 1);
                if query[query_index - 1].base == self.poa_graph.raw_nodes()[topo_indices[graph_index - 1].index()].weight.base {
                    temp_match_score = poa_matrix[temp_match_position][query_index - 1].score + (self.match_score * temp_freq_match as i32) as isize;
                    poa_matrix[graph_index][query_index].back = 'm';
                    //print!("position {} {} mat del in {} {} {} ", graph_index, query_index, temp_match_score - self.match_score as isize, temp_del_score - self.gap_open_score as isize, temp_ins_score - self.gap_open_score as isize);
                    
                }
                else {
                    temp_match_score = poa_matrix[temp_match_position][query_index - 1].score + (self.mismatch_score * temp_freq_query as i32) as isize;
                    poa_matrix[graph_index][query_index].back = 's';
                    //print!("position {} {} mat del in {} {} {} ", graph_index, query_index, temp_match_score - self.mismatch_score as isize, temp_del_score - self.gap_open_score as isize, temp_ins_score - self.gap_open_score as isize);
                }
                // filling out the match matrix
                
                //println!(" real mat del in {} {} {}", temp_match_score , temp_del_score, temp_ins_score);

                if (temp_ins_score >= temp_match_score) && (temp_ins_score >= temp_del_score) {
                    poa_matrix[graph_index][query_index].score = temp_ins_score;
                    poa_matrix[graph_index][query_index].prev = temp_ins_position;
                    poa_matrix[graph_index][query_index].back = 'i';
                    poa_matrix[graph_index][query_index].freq = temp_freq_graph;
                }
                else if temp_del_score >= temp_match_score {
                    poa_matrix[graph_index][query_index].score = temp_del_score;
                    poa_matrix[graph_index][query_index].prev = query_index - 1;
                    poa_matrix[graph_index][query_index].back = 'd';
                    poa_matrix[graph_index][query_index].freq = temp_freq_query;
                }
                else {
                    //print!(" {} == {} ", query[query_index - 1], self.poa_graph.raw_nodes()[topo_indices[graph_index - 1].index()].weight);
                    poa_matrix[graph_index][query_index].score = temp_match_score;
                    poa_matrix[graph_index][query_index].prev = temp_match_position;
                    poa_matrix[graph_index][query_index].freq = temp_freq_match;
                }
            }
        }
         /* 
        println!("printing score matrix");
        print!("{:>3} {:>3} ", 'S', 'S');
        for query_index in 0..query.len() {
            print!("{:>3} ", query[query_index].base as char);
        }
        println!("");
        for graph_index in 0..topo_indices.len() + 1 {
            if graph_index != 0 {
                print!("{:>3} ", self.poa_graph.raw_nodes()[topo_indices[graph_index - 1].index()].weight.base as char);
            }
            else {
                print!("{:>3} ", 'S');
            }
            for query_index in 0..query.len() + 1 {
                print!("{:>3} ", poa_matrix[graph_index][query_index].score);
            }
            println!("");
        }
        println!("printing back matrix");
        print!("{:>3} {:>3} ", 'S', 'S');
        for query_index in 0..query.len() {
            print!("{:>3} ", query[query_index].base as char);
        }
        println!("");
        for graph_index in 0..topo_indices.len() + 1 {
            if graph_index != 0 {
                print!("{:>3} ", self.poa_graph.raw_nodes()[topo_indices[graph_index - 1].index()].weight.base as char);
            }
            else {
                print!("{:>3} ", 'S');
            }
            for query_index in 0..query.len() + 1 {
                print!("{:>3} ", poa_matrix[graph_index][query_index].back);
            }
            println!("");
        }
        */
        // backtrace
        // back tracing using back matrix and filling out align_vec
        let mut align_vec: Vec<(u8, usize, usize)> = vec![];
        let mut i = topo_indices.len();
        let mut j = query.len();
        let mut break_on_next = false;
        loop {
            match poa_matrix[i][j].back {
                'i' => {
                    align_vec.push(('i' as u8, poa_matrix[i][j].prev, poa_matrix[i][j].freq));
                    i = poa_matrix[i][j].prev;
                },
                'm' => {
                    align_vec.push(('m' as u8, poa_matrix[i][j].prev, poa_matrix[i][j].freq));
                    i = poa_matrix[i][j].prev;
                    j = j - 1;
                },
                's' => {
                    align_vec.push(('s' as u8, poa_matrix[i][j].prev, poa_matrix[i][j].freq));
                    i = poa_matrix[i][j].prev;
                    j = j - 1;
                }
                'd' => {
                    align_vec.push(('d' as u8, poa_matrix[i][j].prev, poa_matrix[i][j].freq));
                    j = j - 1;
                },
                _ => (),
            }
            if break_on_next {
                let pos = align_vec.len() - 1;
                align_vec[pos].1 = MAX_USIZE;
                break;
            }
            if i == 0 && j == 0 {
                break_on_next = true;
            }
        }
        for _base in &align_vec {
            //println!("{} {}", base.0 as char, base.1);
        }
        align_vec
    }
    pub fn consensus(&self) -> (Vec<u8>, Vec<usize>) {
        let mut output: Vec<u8> = vec![];
        let mut topopos: Vec<usize> = vec![];
        let mut topo = Topo::new(&self.poa_graph);
        let mut topo_indices = Vec::new();
        let mut max_index = 0;
        let mut max_score = 0.0;
        while let Some(node) = topo.next(&self.poa_graph) {
            topo_indices.push(node);
            if max_index < node.index(){
                max_index = node.index();
            }
        }
        topo_indices.reverse();
        let mut weight_scores: Vec<i32> = vec![0; max_index + 1];
        let mut scores: Vec<f64> = vec![0.0; max_index + 1];
        let mut next_in_path: Vec<usize> = vec![0; max_index + 1];
        for node in topo_indices{
            let mut best_weight_score_edge: (i32, f64, usize) = (-1 , -1.0, MAX_USIZE);
            let mut neighbour_nodes = self.poa_graph.neighbors_directed(node, Outgoing);
            while let Some(neighbour_node) = neighbour_nodes.next() {
                let mut edges = self.poa_graph.edges_connecting(node, neighbour_node);
                let mut weight: i32 = 0;
                while let Some(edge) = edges.next() {
                    weight += edge.weight().clone();
                }
                let weight_score_edge = (weight, scores[neighbour_node.index()], neighbour_node.index());
                if weight_score_edge > best_weight_score_edge{
                    best_weight_score_edge = weight_score_edge;
                }
            }
            if best_weight_score_edge.0 as f64 + best_weight_score_edge.1 > max_score{
                max_score = best_weight_score_edge.0 as f64 + best_weight_score_edge.1;
            }
            scores[node.index()] = best_weight_score_edge.0 as f64 + best_weight_score_edge.1;
            next_in_path[node.index()] = best_weight_score_edge.2;
            weight_scores[node.index()] = best_weight_score_edge.0;
        }
        let mut pos = scores.iter().position(|&r| r == max_score).unwrap();
        let mut consensus_started: bool = false;
        let weight_average = scores[pos] / scores.len() as f64;
        let weight_threshold = weight_average as i32 / 2;
        while pos != MAX_USIZE {
            if consensus_started == false && weight_scores[pos] < weight_threshold {
                pos = next_in_path[pos];
                continue;
            }
            consensus_started = true;
            topopos.push(pos as usize);
            for _ in 0..self.poa_graph.raw_nodes()[pos].weight.frequency {
                output.push(self.poa_graph.raw_nodes()[pos].weight.base);
            }
            pos = next_in_path[pos];
        }
        (output, topopos)
    }
}