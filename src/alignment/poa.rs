
use petgraph::{graph::NodeIndex, visit::Topo, Directed, Graph, Incoming, Outgoing};
use std::cmp;

pub const MIN_ISIZE: isize = -858_993_459;
pub const MAX_USIZE: usize = 858_993_459;

pub struct Poa {
    pub poa_graph: Graph<u8, i32, Directed, usize>,
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
}
#[derive(Clone)]
pub struct PoaMatrixCell {
    match_score: isize,
    insert_score: isize,
    delete_score: isize,
    match_prev: usize,
    insert_prev: usize,
    back: char,
}

impl Poa {
    pub fn initialize (start_query: &Vec<u8>, match_score: i32, mismatch_score: i32, gap_open_score: i32, gap_extend_score: i32) -> Self {
        // allocate space for the graph
        let mut graph: Graph<u8, i32, Directed, usize> = Graph::with_capacity(start_query.len(), start_query.len() - 1);
        // add nodes and edges to the graph
        let mut prev: NodeIndex<usize> = graph.add_node(start_query[0]);
        let mut node: NodeIndex<usize>;
        for base in start_query.iter().skip(1) {
            node = graph.add_node(*base);
            graph.add_edge(prev, node, 1);
            prev = node;
        }
        Poa {poa_graph: graph, match_score: match_score, mismatch_score: mismatch_score, gap_open_score: gap_open_score, gap_extend_score: gap_extend_score}
    }
    pub fn add_to_poa (&mut self, query: &Vec<u8>) {
        // get topological ordering of the graph
        let mut topo = Topo::new(&self.poa_graph);
        let mut topo_indices = Vec::new();
        while let Some(node) = topo.next(&self.poa_graph) {
            topo_indices.push(node);
        }
        let align_vec = self.get_alignment(query);
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
                    // find the edge between alignment.1 and current node
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
                    // make a new node
                    prev_node = current_node;
                    current_node = Some(self.poa_graph.add_node(query[query_index - 1]));
                    match prev_node {
                        Some(x) => {
                            // connect the new node to the current one
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
                    current_node = Some(self.poa_graph.add_node(query[query_index - 1]));
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

    pub fn get_alignment (&self, query: &Vec<u8>) -> Vec<(u8, usize)> {
        // get topological ordering of the graph
        let mut topo = Topo::new(&self.poa_graph);
        let mut topo_indices = Vec::new();
        while let Some(node) = topo.next(&self.poa_graph) {
            topo_indices.push(node);
        }
        // initialize the weight matrix with zeros
        let mut poa_matrix: Vec<Vec<PoaMatrixCell>> = vec![vec![PoaMatrixCell { match_score: (0), insert_score: (0), delete_score: (0), match_prev: (0), insert_prev: (0), back: ('0') }; query.len() + 1]; topo_indices.len() + 1];
        // fill the first row and column
        let mut graph_index = 1;
        let mut query_index = 1;
        // the first column
        for graph_node in &topo_indices {
            if graph_index == 1 {
                let temp_value = self.gap_open_score as isize + self.gap_extend_score as isize;
                poa_matrix[graph_index][0] = PoaMatrixCell { match_score: (temp_value), insert_score: (temp_value), delete_score: (temp_value), match_prev: (0), insert_prev: (0), back: ('i') };
            }
            else if graph_index > 1 {
                let mut max_score = MIN_ISIZE;
                let mut max_position = 0;
                // go through all the incoming nodes to this one and find the one with the max score
                let prevs: Vec<NodeIndex<usize>> = self.poa_graph.neighbors_directed(*graph_node, Incoming).collect();
                for prev_node in &prevs {
                    // find the index of the node in topo list
                    let position = topo_indices.iter().position(|r| r == prev_node).unwrap();
                    // get the score and add gap_extend
                    let temp_score = poa_matrix[position + 1][0].insert_score + self.gap_extend_score as isize;
                    //save the max score and position
                    if temp_score > max_score {
                        max_score = temp_score;
                        max_position = position + 1;
                    }
                }
                poa_matrix[graph_index][0] = PoaMatrixCell { match_score: (max_score), insert_score: (max_score), delete_score: (max_score), match_prev: (max_position), insert_prev: (max_position), back: ('i') };
            }
            graph_index += 1;
        }
        // the first row
        for _ in query {
            if query_index == 1 {
                let temp_value = self.gap_open_score as isize + self.gap_extend_score as isize;
                poa_matrix[0][query_index] = PoaMatrixCell { match_score: (temp_value), insert_score: (temp_value), delete_score: (temp_value), match_prev: (0), insert_prev: (0), back: ('d') };
            }
            else if query_index > 1 {
                let temp_value = poa_matrix[0][query_index - 1].delete_score + self.gap_extend_score as isize;
                poa_matrix[0][query_index] = PoaMatrixCell { match_score: (temp_value), insert_score: (temp_value), delete_score: (temp_value), match_prev: (0), insert_prev: (0), back: ('d') };
            }
            query_index += 1;
        }
        // fill the rest
        for graph_index in 1..topo_indices.len() + 1 {
            for query_index in 1..query.len() + 1 {
                // fill the delete matrix simple
                let temp_del_score = poa_matrix[graph_index][query_index - 1].delete_score + self.gap_extend_score as isize;
                let temp_match_score = poa_matrix[graph_index][query_index - 1].match_score + self.gap_open_score as isize + self.gap_extend_score as isize;
                poa_matrix[graph_index][query_index].delete_score = cmp::max(temp_del_score, temp_match_score);
                // fill the insert matrix complex
                let mut temp_ins_score = MIN_ISIZE;
                let mut temp_ins_position = 0;
                let mut temp_match_score = MIN_ISIZE;
                let mut temp_match_position = 0;
                let prevs: Vec<NodeIndex<usize>> = self.poa_graph.neighbors_directed(topo_indices[graph_index - 1], Incoming).collect();
                for prev_node in &prevs {
                    // find the index of the node in topo list
                    let position = topo_indices.iter().position(|r| r == prev_node).unwrap();
                    // highest insert score and location
                    // get the score and add gap_extend
                    let temp_score = poa_matrix[position + 1][query_index].insert_score + self.gap_extend_score as isize;
                    // save the max score and position
                    if temp_score > temp_ins_score {
                        temp_ins_score = temp_score;
                        temp_ins_position = position + 1;
                    }
                    // highest match score and location
                    let temp_score = poa_matrix[position + 1][query_index].match_score + self.gap_open_score as isize + self.gap_extend_score as isize;
                    // save the max score and position
                    if temp_score > temp_match_score {
                        temp_match_score = temp_score;
                        temp_match_position = position + 1;
                    }
                }
                if graph_index == 1 {
                    poa_matrix[graph_index][query_index].insert_score = poa_matrix[graph_index - 1][query_index].insert_score + self.gap_extend_score as isize;
                    poa_matrix[graph_index][query_index].insert_prev = graph_index - 1;
                }
                else if temp_match_score > temp_ins_score {
                    poa_matrix[graph_index][query_index].insert_score = temp_match_score;
                    poa_matrix[graph_index][query_index].insert_prev = temp_match_position;
                }
                else if temp_ins_score >= temp_match_score {
                    poa_matrix[graph_index][query_index].insert_score = temp_ins_score;
                    poa_matrix[graph_index][query_index].insert_prev = temp_ins_position;
                }
                // fill the match matrix bit more complex ....
                let temp_ins_score = poa_matrix[graph_index][query_index].insert_score;
                // get the i,j from the deletion matrix
                let temp_del_score = poa_matrix[graph_index][query_index].delete_score;
                // match or substitution
                if query[query_index - 1] == self.poa_graph.raw_nodes()[topo_indices[graph_index - 1].index()].weight {
                    temp_match_score = poa_matrix[temp_match_position][query_index - 1].match_score + self.match_score as isize;
                    poa_matrix[graph_index][query_index].back = 'm';
                }
                else {
                    temp_match_score = poa_matrix[temp_match_position][query_index - 1].match_score + self.mismatch_score as isize;
                    poa_matrix[graph_index][query_index].back = 's';
                }
                // filling out the match matrix
                if (temp_match_score >= temp_ins_score) && (temp_match_score >= temp_del_score) {
                    poa_matrix[graph_index][query_index].match_score = temp_match_score;
                    poa_matrix[graph_index][query_index].match_prev = temp_match_position;
                }
                else if temp_ins_score > temp_del_score {
                    poa_matrix[graph_index][query_index].match_score = poa_matrix[graph_index][query_index].insert_score;
                    poa_matrix[graph_index][query_index].match_prev = poa_matrix[graph_index][query_index].insert_prev;
                    poa_matrix[graph_index][query_index].back = 'i';
                }
                else {
                    poa_matrix[graph_index][query_index].match_score = poa_matrix[graph_index][query_index].delete_score;
                    poa_matrix[graph_index][query_index].match_prev = query_index - 1;
                    poa_matrix[graph_index][query_index].back = 'd';
                }
            }
        }
        /* 
        println!("printing ins matrix");
        for graph_index in 0..topo_indices.len() + 1 {
            for query_index in 0..query.len() + 1 {
                print!("{:>3} ", ins_matrix[graph_index][query_index].0);
            }
            println!("");
        }
        println!("printing match matrix");
        for graph_index in 0..topo_indices.len() + 1 {
            for query_index in 0..query.len() + 1 {
                print!("{:>3} ", match_matrix[graph_index][query_index].0);
            }
            println!("");
        }
        println!("printing del matrix");
        for graph_index in 0..topo_indices.len() + 1 {
            for query_index in 0..query.len() + 1 {
                print!("{:>3} ", del_matrix[graph_index][query_index].0);
            }
            println!("");
        }
        for graph_index in 0..topo_indices.len() + 1 {
            for query_index in 0..query.len() + 1 {
                print!("{:>3} ", back_matrix[graph_index][query_index].0);
            }
            println!("");
        }
        */
        // backtrace
        // back tracing using back matrix and filling out align_vec
        let mut align_vec: Vec<(u8, usize)> = vec![];
        let mut i = topo_indices.len();
        let mut j = query.len();
        let mut break_on_next = false;
        loop {
            match poa_matrix[i][j].back {
                'i' => {
                    align_vec.push(('i' as u8, poa_matrix[i][j].insert_prev));
                    i = poa_matrix[i][j].insert_prev;
                },
                'm' => {
                    align_vec.push(('m' as u8, poa_matrix[i][j].match_prev));
                    i = poa_matrix[i][j].match_prev;
                    j = j - 1;
                },
                's' => {
                    align_vec.push(('s' as u8, poa_matrix[i][j].match_prev));
                    i = poa_matrix[i][j].match_prev;
                    j = j - 1;
                }
                'd' => {
                    align_vec.push(('d' as u8, i));
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
        for base in &align_vec {
            println!("{} {}", base.0 as char, base.1);
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
        //define score and nextinpath vectors with capacity of num nodes.
        let mut weight_scores: Vec<i32> = vec![0; max_index + 1];
        let mut scores: Vec<f64> = vec![0.0; max_index + 1];
        let mut next_in_path: Vec<usize> = vec![0; max_index + 1];
        //iterate thorugh the nodes in revere
        for node in topo_indices{
            //print!("\nstart node: {:?}", self.graph.raw_nodes()[node.index()].weight);
            let mut best_weight_score_edge: (i32, f64, usize) = (-1 , -1.0, 123456789);
            //let mut outEdges = self.graph.neighbors_directed(node, Outgoing).detach();
            let mut neighbour_nodes = self.poa_graph.neighbors_directed(node, Outgoing);
            while let Some(neighbour_node) = neighbour_nodes.next() {
                //print!(" end node: {:?}", self.graph.raw_nodes()[neighbour_node.index()].weight);
                let mut edges = self.poa_graph.edges_connecting(node, neighbour_node);
                let mut weight: i32 = 0;
                while let Some(edge) = edges.next() {
                    weight += edge.weight().clone();
                    //print!(" Edge of weight {}", weight);
                }
                let weight_score_edge = (weight, scores[neighbour_node.index()], neighbour_node.index());
                if weight_score_edge > best_weight_score_edge{
                    best_weight_score_edge = weight_score_edge;
                }
            }
            //save score and traceback
            if best_weight_score_edge.0 as f64 + best_weight_score_edge.1 > max_score{
                max_score = best_weight_score_edge.0 as f64 + best_weight_score_edge.1;
            }
            scores[node.index()] = best_weight_score_edge.0 as f64 + best_weight_score_edge.1;
            next_in_path[node.index()] = best_weight_score_edge.2;
            weight_scores[node.index()] = best_weight_score_edge.0;
        }
        let mut pos = scores.iter().position(|&r| r == max_score).unwrap();
        //calculate the start weight score
        let mut consensus_started: bool = false;
        let weight_average = scores[pos] / scores.len() as f64;
        let weight_threshold = weight_average as i32 / 2;
        while pos != 123456789 {
            //continue if starting weight score is too low
            if consensus_started == false && weight_scores[pos] <= weight_threshold {
                pos = next_in_path[pos];
                continue;
            }
            consensus_started = true;
            topopos.push(pos as usize);
            output.push(self.poa_graph.raw_nodes()[pos].weight);
            pos = next_in_path[pos];
        }
        (output, topopos)
    }
}