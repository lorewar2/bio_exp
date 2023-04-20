
use petgraph::{Direction::Outgoing, graph::NodeIndex, visit::Topo, Directed, Graph, Incoming};
use std::cmp;

pub const MIN_ISIZE: isize = -858_993_459;
pub struct Poa {
    poa_graph: Graph<u8, i32, Directed, usize>,
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
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
    pub fn add_to_poa (query: &Vec<u8>) {
        
    }

    pub fn get_alignment (self, query: &Vec<u8>) {
        // get topological ordering of the graph
        let mut topo = Topo::new(&self.poa_graph);
        let mut topo_indices = Vec::new();
        while let Some(node) = topo.next(&self.poa_graph) {
            topo_indices.push(node);
        }
        // make three matrices
        // initialize the weight matrix with zeros
        let mut match_matrix: Vec<Vec<(isize, (usize, usize))>> = vec![vec![(0, (0, 0)); query.len() + 1]; self.poa_graph.node_count() + 1]; // match or mismatch diagonal edges  // score and the aligned one location
        let mut del_matrix: Vec<Vec<(isize, (usize, usize))>> = vec![vec![(0, (0, 0)); query.len() + 1]; self.poa_graph.node_count() + 1];  // x deletions right direction edges // score and the aligned one location
        let mut ins_matrix: Vec<Vec<(isize, (usize, usize))>> = vec![vec![(0, (0, 0)); query.len() + 1]; self.poa_graph.node_count() + 1]; // x insertion down direction edges // score and the aligned one location
        let mut back_matrix: Vec<Vec<(char, (usize, usize))>> = vec![vec![('m', (0, 0)); query.len() + 1]; self.poa_graph.node_count() + 1];
        // fill the first row and column
        let mut graph_index = 1;
        let mut query_index = 1;
        // the first column
        for graph_node in &topo_indices {
            if graph_index == 1 {
                back_matrix[graph_index][0] = ('i', (0, 0));
                ins_matrix[graph_index][0] = (self.gap_open_score as isize + self.gap_extend_score as isize, (0, 0));
                del_matrix[graph_index][0] = (self.gap_open_score as isize + self.gap_extend_score as isize, (0, 0));
                match_matrix[graph_index][0] = (self.gap_open_score as isize + self.gap_extend_score as isize, (0, 0));
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
                    let temp_score = ins_matrix[position + 1][0].0 + self.gap_extend_score as isize;
                    //save the max score and position
                    if temp_score > max_score {
                        max_score = temp_score;
                        max_position = position + 1;
                    }
                }
                back_matrix[graph_index][0] = ('i', (max_position, 0)); // max position is the matrix position not sequence position (seq_pos = max_pos - 1)
                ins_matrix[graph_index][0] = (max_score, (max_position, 0));
                del_matrix[graph_index][0] = (max_score, (max_position, 0));
                match_matrix[graph_index][0] = (max_score, (max_position, 0));
            }
            graph_index += 1;
        }
        // the first row
        for _ in query {
            if query_index == 1 {
                back_matrix[0][query_index] = ('d', (0, 0));
                ins_matrix[0][query_index] = (self.gap_open_score as isize + self.gap_extend_score as isize, (0, 0));
                del_matrix[0][query_index] = (self.gap_open_score as isize + self.gap_extend_score as isize, (0, 0));
                match_matrix[0][query_index] = (self.gap_open_score as isize + self.gap_extend_score as isize, (0, 0));
            }
            else if query_index > 1 {
                back_matrix[0][query_index] = ('d', (0, query_index - 1));
                del_matrix[0][query_index] = (del_matrix[0][query_index - 1].0 + self.gap_extend_score as isize, (0, query_index - 1));
                ins_matrix[0][query_index] = (del_matrix[0][query_index - 1].0 + self.gap_extend_score as isize, (0, query_index - 1));
                match_matrix[0][query_index] = (del_matrix[0][query_index - 1].0 + self.gap_extend_score as isize, (0, query_index - 1));
            }
            query_index += 1;
        }
        // fill the rest
        for graph_index in 1..topo_indices.len() + 1 {
            for query_index in 1..query.len() + 1 {
                // fill the delete matrix simple
                let temp_del_score = del_matrix[graph_index][query_index - 1].0 + self.gap_extend_score as isize;
                let temp_match_score = match_matrix[graph_index][query_index - 1].0 + self.gap_open_score as isize + self.gap_extend_score as isize;
                del_matrix[graph_index][query_index] = (cmp::max(temp_del_score, temp_match_score), (graph_index, query_index - 1));
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
                    let temp_score = ins_matrix[position + 1][query_index].0 + self.gap_extend_score as isize;
                    // save the max score and position
                    if temp_score > temp_ins_score {
                        temp_ins_score = temp_score;
                        temp_ins_position = position + 1;
                    }
                    // highest match score and location
                    let temp_score = match_matrix[position + 1][query_index].0 + self.gap_open_score as isize + self.gap_extend_score as isize;
                    // save the max score and position
                    if temp_score > temp_match_score {
                        temp_match_score = temp_score;
                        temp_match_position = position + 1;
                    }
                }
                if graph_index == 1 {
                    ins_matrix[graph_index][query_index] = (ins_matrix[graph_index - 1][query_index].0 + self.gap_extend_score as isize, (graph_index - 1, query_index));
                }
                else if temp_match_score > temp_ins_score {
                    ins_matrix[graph_index][query_index] = (temp_match_score, (temp_match_position, query_index));
                }
                else if temp_ins_score >= temp_match_score {
                    ins_matrix[graph_index][query_index] = (temp_ins_score, (temp_ins_position, query_index));
                }
                // fill the match matrix bit more complex ....
                let temp_ins_score = ins_matrix[graph_index][query_index].0;
                // get the i,j from the deletion matrix
                let temp_del_score = del_matrix[graph_index][query_index].0;
                // match or substitution
                if query[query_index - 1] == self.poa_graph.raw_nodes()[topo_indices[graph_index - 1].index()].weight {
                    temp_match_score = match_matrix[temp_match_position][query_index - 1].0 + self.match_score as isize;
                    back_matrix[graph_index][query_index] = ('m', (temp_match_position, query_index - 1));
                }
                else {
                    temp_match_score = match_matrix[temp_match_position][query_index - 1].0 + self.mismatch_score as isize;
                    back_matrix[graph_index][query_index] = ('s', (temp_match_position, query_index - 1));
                }
                // filling out the match matrix
                if (temp_match_score >= temp_ins_score) && (temp_match_score >= temp_del_score) {
                    match_matrix[graph_index][query_index] = (temp_match_score, (temp_match_position, query_index));
                }
                else if temp_ins_score > temp_del_score {
                    match_matrix[graph_index][query_index] = ins_matrix[graph_index][query_index];
                    back_matrix[graph_index][query_index] = ('i', ins_matrix[graph_index][query_index].1);
                }
                else {
                    match_matrix[graph_index][query_index] = del_matrix[graph_index][query_index];
                    back_matrix[graph_index][query_index] = ('d', del_matrix[graph_index][query_index].1);
                }
            }
        }
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
        // backtrace
        // back tracing using back matrix and filling out align_vec
    let mut align_vec: Vec<(u8, usize, usize)> = vec![];
    let mut i = topo_indices.len();
    let mut j = query.len();
    let score = match_matrix[i][j].0;
    loop {
        match back_matrix[i][j].0 {
            'i' => {
                align_vec.push(('i' as u8, back_matrix[i][j].1.0, back_matrix[i][j].1.1));
                i = back_matrix[i][j].1.0;
                
            },
            'm' => {
                align_vec.push(('m' as u8, back_matrix[i][j].1.0, back_matrix[i][j].1.1));
                i = back_matrix[i][j].1.0;
                j = j - 1;
                
            },
            's' => {
                align_vec.push(('s' as u8, back_matrix[i][j].1.0, back_matrix[i][j].1.1));
                i = back_matrix[i][j].1.0;
                j = j - 1;
                
            }
            'd' => {
                align_vec.push(('d' as u8, back_matrix[i][j].1.0, back_matrix[i][j].1.1));
                j = j - 1;
               
            },
            _ => (),
        }
        if i == 0 && j == 0 {
            break;
        }
    }
    println!("{:?}", align_vec);
    }
}