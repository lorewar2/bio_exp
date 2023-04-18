
use petgraph::{Direction::Outgoing, graph::NodeIndex, visit::Topo, Directed, Graph, Incoming};
use std::cmp;

pub const MIN_ISIZE: isize = -858_993_459;
pub struct poa {
    poa_graph: Graph<u8, i32, Directed, usize>,
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
}

impl poa {
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
        poa {poa_graph: graph, match_score: match_score, mismatch_score: mismatch_score, gap_open_score: gap_open_score, gap_extend_score: gap_extend_score}
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
        topo = Topo::new(&self.poa_graph);
        // make three matrices
        // initialize the weight matrix with zeros
        let mut match_matrix: Vec<Vec<(isize, (usize, usize))>> = vec![vec![(0, (0, 0)); query.len() + 1]; self.poa_graph.node_count() + 1]; // match or mismatch diagonal edges  // score and the aligned one location
        let mut del_matrix: Vec<Vec<(isize, (usize, usize))>> = vec![vec![(0, (0, 0)); query.len() + 1]; self.poa_graph.node_count() + 1];  // x deletions right direction edges // score and the aligned one location
        let mut ins_matrix: Vec<Vec<(isize, (usize, usize))>> = vec![vec![(0, (0, 0)); query.len() + 1]; self.poa_graph.node_count() + 1]; // x insertion down direction edges // score and the aligned one location
        let mut back_matrix: Vec<Vec<(char, (usize, usize))>> = vec![vec![('m', (0, 0)); query.len() + 1]; self.poa_graph.node_count() + 1];
        // fill the first row and column
        let mut graph_index = 0;
        let mut query_index = 0;
        // the first column
        while let Some(graph_node) = topo.next(&self.poa_graph) {
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
                let prevs: Vec<NodeIndex<usize>> = self.poa_graph.neighbors_directed(graph_node, Incoming).collect();
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
                
            }
            graph_index += 1;
        }
        // the first row

        // fill the rest

        // backtrace
    }
}