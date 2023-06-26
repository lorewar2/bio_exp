use std::cmp;

use petgraph::{graph::{NodeIndex}, visit::Topo, Directed, Graph, Incoming, Outgoing, Direction};
use libm::exp;
use logaddexp::LogAddExp;
use statrs::function::factorial::binomial;

const ERROR_PROBABILITY: f64 = 0.80;
const PRINT_ALL: bool = false;
const USEPACBIODATA: bool = true;
const NUM_OF_ITER_FOR_PARALLEL: usize = 10;

pub fn get_consensus_quality_scores(seq_num: usize, consensus: &Vec<u8>, topology: &Vec<usize>, graph: &Graph<u8, i32, Directed, usize>) -> (Vec<f64>, Vec<Vec<usize>>) {
    let mut quality_scores: Vec<f64> = vec![];
    let mut base_count_vec: Vec<Vec<usize>> = vec![];
    //run all the consensus through get indices
    for i in 0..consensus.len() {
        println!("{} / {}", i, consensus.len());
        // skip the indices which are in the passed consensus
        let skip_nodes: Vec<usize> = topology[0 .. i + 1].to_vec();
        // new method using topology cut
        let mut target_node_parent = None;
        let mut target_node_child = None;
        if i != 0{
            target_node_parent = Some(topology[i - 1]);
        }
        if i != consensus.len() - 1 {
            target_node_child = Some(topology[i + 1]);
        }
        let (parallel_nodes, parallel_num_incoming_seq, _) = get_parallel_nodes_with_topology_cut (skip_nodes, seq_num,  topology[i], target_node_parent, target_node_child, graph);
        let (temp_quality_score, _, temp_base_counts, _) = base_quality_score_calculation (seq_num, parallel_nodes, parallel_num_incoming_seq, consensus[i], graph);
        quality_scores.push(temp_quality_score);
        base_count_vec.push(temp_base_counts);
    }
    (quality_scores, base_count_vec)
}

pub fn get_parallel_nodes_with_topology_cut (skip_nodes: Vec<usize>, total_seq: usize, target_node: usize, target_node_parent: Option<usize>, target_node_child: Option<usize>, graph: &Graph<u8, i32, Directed, usize>) -> (Vec<usize>, Vec<usize>, Vec<String>) {
    // vector initialization
    let mut debug_strings: Vec<String> = vec![];
    let mut topology = Topo::new(graph);
    let mut topologically_ordered_nodes = Vec::new();
    let mut parallel_nodes: Vec<usize> = vec![];
    let mut parallel_node_parents: Vec<usize> = vec![];
    let mut parallel_num_incoming_seq: Vec<usize> = vec![];
    let mut direction: Option<Direction> = None;
    let temp_string;
    //print stuff
    if PRINT_ALL {
        println!("NODE CHECKING FOR PARALLEL: {}, base {}", target_node, graph.raw_nodes()[target_node].weight as char);
    }

    // make a topologically ordered list
    while let Some(node) = topology.next(graph) {
        topologically_ordered_nodes.push(node.index());
    }
    // find the position of the target node, its child and parent in topology list
    let target_node_topological_position = topologically_ordered_nodes.iter().position(|&r| r == target_node).unwrap();
    let target_child_topological_position = match target_node_child { 
        Some(child_node) => {topologically_ordered_nodes.iter().position(|&r| r == child_node).unwrap()},
        None => {direction = Some(Incoming); target_node_topological_position}
    };
    let target_parent_topological_position = match target_node_parent { 
        Some(parent_node) => {topologically_ordered_nodes.iter().position(|&r| r == parent_node).unwrap()},
        None => {direction = Some(Outgoing); target_node_topological_position}
    };
    // choose a direction with the least amount of intermediate nodes
    if (direction == None) && (topologically_ordered_nodes[target_parent_topological_position..target_node_topological_position].len() > topologically_ordered_nodes[target_node_topological_position..target_child_topological_position].len()) {
        direction = Some(Outgoing);
    }
    else if direction == None {
        direction = Some(Incoming);
    }
    match direction {
        Some(x) => {
            if x == Incoming {
                temp_string = format!("Going backwards");
                if PRINT_ALL {
                    println!("{}", temp_string);
                }
                debug_strings.push(temp_string.clone());
            }
            else {
                temp_string = format!("Going forward");
                if PRINT_ALL {
                    println!("{}", temp_string);
                }
                debug_strings.push(temp_string.clone());
            }
        }
        None => {}
    }
    // check if the target node corrosponds with all the sequences
    let num_seq_through_target_base = find_the_seq_passing_through (target_node, graph);

    if num_seq_through_target_base == total_seq {
        parallel_nodes.push(target_node);
        if USEPACBIODATA {
            parallel_num_incoming_seq.push(num_seq_through_target_base - 1);
        }
        parallel_num_incoming_seq.push(num_seq_through_target_base);
        return (parallel_nodes, parallel_num_incoming_seq, debug_strings);
    }
    // go back skip_count and go forward skip_count + 3 and check if parent and child are before and after target_node_position,
    // iterate skip_count until all sequences are found, break on 5
    let mut seq_found_so_far = num_seq_through_target_base;
    let mut bubble_size = 1;
    while (seq_found_so_far < total_seq)  && (bubble_size < NUM_OF_ITER_FOR_PARALLEL) {
        let temp_debug_strings;
        (parallel_nodes, parallel_node_parents, parallel_num_incoming_seq, seq_found_so_far, temp_debug_strings) = move_in_direction_and_find_crossing_nodes (&skip_nodes, total_seq, direction.unwrap(), parallel_nodes, parallel_node_parents, parallel_num_incoming_seq, seq_found_so_far, target_node, bubble_size, &topologically_ordered_nodes, target_node_topological_position, graph);
        debug_strings = [debug_strings, temp_debug_strings].concat();
        bubble_size += 1;
    }
    if USEPACBIODATA {
        parallel_num_incoming_seq.push(num_seq_through_target_base - 1);
    }
    parallel_num_incoming_seq.push(num_seq_through_target_base);
    parallel_nodes.push(target_node);
    (parallel_nodes, parallel_num_incoming_seq, debug_strings)
}

fn find_the_seq_passing_through (target_node: usize, graph: &Graph<u8, i32, Directed, usize>) -> usize {
    let node_index = NodeIndex::new(target_node);
    //edges directed toward the base
    let incoming_nodes: Vec<NodeIndex<usize>> = graph.neighbors_directed(node_index, Incoming).collect();
    let mut incoming_weight = 0;
    for incoming_node in incoming_nodes {
        let mut edges = graph.edges_connecting(incoming_node, node_index);
        while let Some(edge) = edges.next() {
            incoming_weight += edge.weight().clone();
        }
    }
    //edges directed from the base
    let outgoing_nodes: Vec<NodeIndex<usize>> = graph.neighbors_directed(node_index, Outgoing).collect();
    let mut outgoing_weight = 0;
    for outgoing_node in outgoing_nodes {
        let mut edges = graph.edges_connecting(node_index, outgoing_node);
        while let Some(edge) = edges.next() {
            outgoing_weight += edge.weight().clone();
        }
    }
    cmp::max(outgoing_weight, incoming_weight) as usize
}

fn move_in_direction_and_find_crossing_nodes (skip_nodes: &Vec<usize>, total_seq: usize, direction: Direction, mut parallel_nodes: Vec<usize>, mut parallel_node_parents: Vec<usize>, mut parallel_num_incoming_seq: Vec<usize>, mut seq_found_so_far: usize, focus_node: usize, bubble_size: usize, topologically_ordered_nodes: &Vec<usize>, target_node_position: usize, graph: &Graph<u8, i32, Directed, usize>) -> (Vec<usize>, Vec<usize>, Vec<usize>, usize, Vec<String>) {
    let mut debug_strings: Vec<String> = vec![];
    let mut temp_string: String;
    // get a list of x back_iterations back nodes
    let back_nodes_list = get_xiterations_direction_nodes(direction, bubble_size, vec![], focus_node, graph);
    // get a list of all forward nodes 0..(back_iterations + 3) for all the back_nodes
    let mut edge_nodes_list: Vec<usize> = vec![];
    for back_node in &back_nodes_list {
        let temp_forward_list = get_direction_nodes (direction.opposite(), bubble_size + 3, vec![], *back_node, graph);
        for temp_forward_node in &temp_forward_list {
            if !edge_nodes_list.contains(temp_forward_node) {
                edge_nodes_list.push(*temp_forward_node);
            }
        }
    }
    temp_string = format!("Iteration: {} BackNodes: {:?} CheckNodes: {:?}", bubble_size, back_nodes_list, edge_nodes_list);
    if PRINT_ALL {
        println!("{}", temp_string);
    }
    debug_strings.push(temp_string.clone());
    
    // get the two slices of topologically_ordered_list back front
    let mut slice: Vec<Vec<usize>> = [topologically_ordered_nodes[0..target_node_position].to_vec(), topologically_ordered_nodes[target_node_position + 1..topologically_ordered_nodes.len()].to_vec()].to_vec();
    // for debugging
    if slice[0].len() > 10 {
        temp_string = format!("Back slice {:?}", slice[0][(slice[0].len() - 10)..slice[0].len()].to_vec());
        if PRINT_ALL {
            println!("{}", temp_string);
        }
        debug_strings.push(temp_string.clone());
    }
    else {
        temp_string = format!("Back slice {:?}", slice[0][0..slice[0].len()].to_vec());
        if PRINT_ALL {
            println!("{}", temp_string);
        }
        debug_strings.push(temp_string.clone());
    }
    if slice[1].len() > 10 {
        temp_string = format!("Front slice {:?}", slice[1][0..10].to_vec());
        if PRINT_ALL {
            println!("{}", temp_string);
        }
        debug_strings.push(temp_string.clone());
    }
    else {
        temp_string = format!("Front slice {:?}", slice[1][0..slice[1].len()].to_vec());
        if PRINT_ALL {
            println!("{}", temp_string);
        }
        debug_strings.push(temp_string.clone());
    }

    if direction == Outgoing {
        slice.reverse();
    }
    //iterate through edge nodes obtained
    for edge_node in &edge_nodes_list {
        // get the parents of the edge node
        let edge_node_parents = get_direction_nodes (direction, 1, vec![], *edge_node, graph);
        'parent_loop: for edge_node_parent in &edge_node_parents {
            // if the parent is in back section and node is in front section add to parallel nodes or if both parent and target is in intermediate add to parallel loop
            if slice[0].contains(edge_node_parent) && slice[1].contains(edge_node) && (*edge_node_parent != focus_node) {
                // edge node parent check
                if parallel_nodes.contains(edge_node) && parallel_node_parents.contains(edge_node_parent) {
                    // go through the parallel nodes and if there is a match check if the same parent and continue if so
                    for index in 0..parallel_nodes.len() {
                        if (parallel_nodes[index] == *edge_node) && (parallel_node_parents[index] == *edge_node_parent) {
                            continue 'parent_loop;
                        }
                    }
                }
                // target node front of parallel node check
                if direction == Incoming {
                    if get_direction_nodes(Outgoing, 4, vec![],  *edge_node, graph).contains(&focus_node) {
                        continue;
                    }
                    if skip_nodes.contains(edge_node) {
                        continue;
                    }
                }
                else {
                    if get_direction_nodes(Outgoing, 4, vec![], focus_node, graph).contains(&edge_node) {
                        continue;
                    }
                    if skip_nodes.contains(edge_node_parent) {
                        continue;
                    }
                }
                // all found 
                if seq_found_so_far >= total_seq {
                    break;
                }
                parallel_nodes.push(*edge_node);
                parallel_node_parents.push(*edge_node_parent);
                temp_string = format!("success node {} parent/child {}\n", *edge_node, *edge_node_parent);
                if PRINT_ALL {
                    println!("{}", temp_string);
                }
                debug_strings.push(temp_string.clone());
                // get the edge weight and add to seq_found_so_far
                let mut incoming_weight = 0;
                if direction == Incoming {
                    let mut edges = graph.edges_connecting(NodeIndex::new(*edge_node_parent), NodeIndex::new(*edge_node));
                    while let Some(edge) = edges.next() {
                        incoming_weight += edge.weight().clone();
                    }
                }
                else {
                    let mut edges = graph.edges_connecting(NodeIndex::new(*edge_node), NodeIndex::new(*edge_node_parent));
                    while let Some(edge) = edges.next() {
                        incoming_weight += edge.weight().clone();
                    }
                }
                parallel_num_incoming_seq.push(incoming_weight as usize);
                seq_found_so_far += incoming_weight as usize;
            }
            
        }
    }
    (parallel_nodes, parallel_node_parents, parallel_num_incoming_seq, seq_found_so_far, debug_strings)
}

fn get_direction_nodes (direction: Direction, iteration: usize, mut direction_node_list: Vec<usize>, focus_node: usize, graph: &Graph<u8, i32, Directed, usize>) -> Vec<usize> {
    //forward outgoing
    //backward incoming
    if iteration <= 0 {
        return direction_node_list;
    }
    //get the back nodes of the target
    let mut direction_neighbours = graph.neighbors_directed(NodeIndex::new(focus_node), direction);
    //iterate through the neighbours
    while let Some(direction_neighbour) = direction_neighbours.next() {
        if !direction_node_list.contains(&direction_neighbour.index()){
            direction_node_list.push(direction_neighbour.index());
            direction_node_list = get_direction_nodes (direction, iteration - 1, direction_node_list, direction_neighbour.index(), graph);
        }
    }
    direction_node_list
}

fn get_xiterations_direction_nodes (direction: Direction ,iteration: usize, mut direction_node_list: Vec<usize>, focus_node: usize, graph: &Graph<u8, i32, Directed, usize>) -> Vec<usize> {
    if iteration <= 0 {
        return direction_node_list;
    }
    //get the back nodes of the target
    let mut direction_neighbours = graph.neighbors_directed(NodeIndex::new(focus_node), direction);
    //iterate through the neighbours
    while let Some(direction_neighbour) = direction_neighbours.next() {
        if iteration == 1 {
            if !direction_node_list.contains(&direction_neighbour.index()){
                direction_node_list.push(direction_neighbour.index());
            }
        }
        direction_node_list = get_xiterations_direction_nodes (direction, iteration - 1, direction_node_list, direction_neighbour.index(), graph);
    }
    direction_node_list
}

pub fn base_quality_score_calculation (mut total_seq: usize, indices_of_parallel_nodes: Vec<usize>, seq_through_parallel_nodes: Vec<usize>, base: u8, graph: &Graph<u8, i32, Directed, usize>) -> (f64, bool, Vec<usize>, Vec<String>) {
    if USEPACBIODATA {
        total_seq -= 1;
    }
    //variable initialization
    let mut debug_strings: Vec<String> = Vec::new();
    let mut count_mismatch: bool = false;
    let error_score: f64;
    let quality_score;
    let base_counts: Vec<usize>;

    let ln_prob_base_a = 0.25_f64.ln();
    let ln_prob_base_c = 0.25_f64.ln();
    let ln_prob_base_g = 0.25_f64.ln();
    let ln_prob_base_t = 0.25_f64.ln();
    
    let mut base_a_count = 0;
    let mut base_c_count = 0;
    let mut base_g_count = 0;
    let mut base_t_count = 0;
    //find out how many sequences run through each base
    //match the indices to the base and ++
    for index in 0..indices_of_parallel_nodes.len() {
        match graph.raw_nodes()[indices_of_parallel_nodes[index]].weight {
            65 => {
                base_a_count += seq_through_parallel_nodes[index];
            },
            67 => {
                base_c_count += seq_through_parallel_nodes[index];
            },
            71 => {
                base_g_count += seq_through_parallel_nodes[index];
            },
            84 => {
                base_t_count += seq_through_parallel_nodes[index];
            },
            _ => {
                //nothing
                },
        }
    }
    // save the base counts for debug
    base_counts = [base_a_count, base_c_count, base_g_count, base_t_count].to_vec();
    if (base_a_count + base_c_count + base_g_count + base_t_count) != (total_seq) {
        count_mismatch = true;
    }
    match count_mismatch {
        true => {
            let temp_string = format!("base counts A:{} C:{} G:{} T:{} MISMATCHHHH!!!!!!!!!!!!!!!!!!!!! \n", base_a_count, base_c_count, base_g_count, base_t_count);
            //println!("{}", temp_string);
            debug_strings.push(temp_string.clone());
        },
        false => {
            let temp_string = format!("base counts A:{} C:{} G:{} T:{}\n", base_a_count, base_c_count, base_g_count, base_t_count);
            //println!("{}", temp_string);
            debug_strings.push(temp_string.clone());
        }
    }
    
    //calculate all the probablilities
    let ln_prob_data_given_a = calculate_binomial(total_seq, base_a_count, ERROR_PROBABILITY).ln();
    let ln_prob_data_given_c = calculate_binomial(total_seq, base_c_count, ERROR_PROBABILITY).ln();
    let ln_prob_data_given_g = calculate_binomial(total_seq, base_g_count, ERROR_PROBABILITY).ln();
    let ln_prob_data_given_t = calculate_binomial(total_seq, base_t_count, ERROR_PROBABILITY).ln();
    //println!("D|A:{} D|C:{} D|G:{} D|T:{}", ln_prob_data_given_a, ln_prob_data_given_c, ln_prob_data_given_g, ln_prob_data_given_t);
    //get the error score *changed to log space*
    let mut ln_sum_of_probablities = ln_prob_data_given_a + ln_prob_base_a;
    ln_sum_of_probablities = ln_sum_of_probablities.ln_add_exp(ln_prob_data_given_c + ln_prob_base_c);
    ln_sum_of_probablities = ln_sum_of_probablities.ln_add_exp(ln_prob_data_given_g + ln_prob_base_g);
    ln_sum_of_probablities = ln_sum_of_probablities.ln_add_exp(ln_prob_data_given_t + ln_prob_base_t);

    error_score =  1.0 - match base {
        65 => {
            //println!("Focus base : A" );
            exp(ln_prob_data_given_a + ln_prob_base_a - ln_sum_of_probablities)
        },
        67 => {
            //println!("Focus base : C" );
            exp(ln_prob_data_given_c + ln_prob_base_c - ln_sum_of_probablities)
        },
        71 => {
            //println!("Focus base : G" );
            exp(ln_prob_data_given_g + ln_prob_base_g - ln_sum_of_probablities)
        },
        84 => {
            //println!("Focus base : T" );
            exp(ln_prob_data_given_t + ln_prob_base_t - ln_sum_of_probablities)
        },
        _ => {0.0},
    };
    quality_score = (-10.00) * error_score.log10();
    //println!("quality score: {}", quality_score);
    //println!("");
    (quality_score, count_mismatch, base_counts, debug_strings)
}

fn calculate_binomial (n: usize, k: usize, prob: f64) -> f64 {
    let binomial_coeff = binomial(n as u64, k as u64);
    let success: f64 = prob.powf(k as f64);
    let failure: f64 = (1.00 - prob).powf((n - k) as f64);
    binomial_coeff * success * failure
}