use crate::alignment::pairwise::pairwise;
use crate::alignment::poabandedsmarter::Aligner;
use crate::generator::simple::get_random_sequences_from_generator;
use crate::alignment::poahomopolymer::Poa;
use crate::quality::topology_cut::base_quality_score_calculation;
use crate::quality::topology_cut::get_parallel_nodes_with_topology_cut;
use crate::quality::topology_cut::get_consensus_quality_scores;
use petgraph::{Graph, Directed, graph::NodeIndex};
use petgraph::dot::Dot;
use rust_htslib::bam::{Read as BamRead, IndexedReader as BamIndexedReader};
use rust_htslib::bcf::{Reader, Read as BcfRead};
use rust_htslib::faidx;
use std::{fs::OpenOptions, io::{prelude::*}};
use std::io::BufReader;
use std::fs::File;
use std::io::SeekFrom;
use std::time::Instant;
use std::fs::read_dir;
use std::fs::create_dir_all;
use std::{io};
use std::fs::read_to_string;

const SEED: u64 = 2;
const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = 0;
const MATCH: i32 = 2;
const MISMATCH: i32 = -2;
const RANDOM_SEQUENCE_LENGTH: usize = 20000;
const NUMBER_OF_RANDOM_SEQUENCES: usize = 5;
const THREE_BASE_CONTEXT_READ_LENGTH: usize = 2;
const NUM_OF_ITER_FOR_ZOOMED_GRAPHS: usize = 4;
const DATA_PATH: &str = "/data1/hifi_consensus/try2/";
const READ_BAM_PATH: &str = "/data1/hifi_consensus/try2/merged.bam";
const INTERMEDIATE_PATH: &str = "result/intermediate";
const CONFIDENT_PATH: &str = "/data1/GiaB_benchmark/HG001_GRCh38_1_22_v4.2.1_benchmark.bed";
const BAND_SIZE: i32 = 2000;
const MAX_NODES_IN_POA: usize = 50000;

pub fn test_graphs() {
    let seqvec = get_random_sequences_from_generator(RANDOM_SEQUENCE_LENGTH, NUMBER_OF_RANDOM_SEQUENCES, SEED);
    println!("Processing seq 1");
    let mut aligner = Aligner::new(MATCH, MISMATCH, GAP_OPEN, &seqvec[0].as_bytes().to_vec(), 200);
    let mut index = 0;
    for seq in &seqvec {
        if index != 0 {
            println!("Processing seq {}", index + 1);
            aligner.global(&seq.as_bytes().to_vec()).add_to_graph();
        }
        index += 1;
    }
    let test_graph = aligner.graph();
    save_the_graph(test_graph, "test.txt".to_string());
    load_the_graph("test.txt".to_string());
}

pub fn save_the_graph (graph: &Graph<u8, i32, Directed, usize>, file_name: String) {
    // check if file is available
    if check_file_availability(&file_name, INTERMEDIATE_PATH) == false {
        let write_string = format!("{}\n{:?}", graph.node_count(), Dot::new(&graph.map(|_, n| (*n) as char, |_, e| *e)));
        let write_path = format!("{}/{}", INTERMEDIATE_PATH, file_name);
        write_string_to_file(&write_path, &write_string);
    }
}

pub fn load_the_graph (file_name: String) -> Graph<u8, i32, Directed, usize> {
    let mut node_edge_list: Vec<(char, Vec<(usize, usize)>)> = vec![];
    let mut node_capacity = 0;
    let mut edge_capacity = 0;
    // check if available, populate node_edge_list from file
    if check_file_availability(&file_name, INTERMEDIATE_PATH) == true {
        // read the file
        let mut index = 0;
        let read_path = format!("{}/{}", INTERMEDIATE_PATH, file_name);
        for line in read_to_string(read_path).unwrap().lines() {
            //println!("{}", line);
            let line_parts: Vec<&str> = line.split(" ").collect();
            // create the node edge list with capacity
            if index == 0 {
                if line_parts.len() == 1 {
                    node_capacity = line_parts[0].parse::<usize>().unwrap();
                    node_edge_list = vec![('X', vec![]); node_capacity];
                }
            }
            // put the values in node edge list
            else {
                if line_parts.len() == 10 {
                    let start_node = line_parts[4].parse::<usize>().unwrap();
                    let chars: Vec<char> = line_parts[8].chars().collect();
                    node_edge_list[start_node].0 = chars[2];
                }
                if line_parts.len() == 12 {
                    let start_node = line_parts[4].parse::<usize>().unwrap();
                    let end_node = line_parts[6].parse::<usize>().unwrap();
                    let chars: Vec<char> = line_parts[10].chars().collect();
                    let edge_weight = chars[1].to_string().parse::<usize>().unwrap();
                    node_edge_list[start_node].1.push((end_node, edge_weight));
                    edge_capacity += 1;
                }

            }
            index += 1;
        }
    }
    // make the graph from the vector
    let mut graph: Graph<u8, i32, Directed, usize> = Graph::with_capacity(node_capacity, edge_capacity);
    
    // add the nodes first
    for node_edge in &node_edge_list {
        graph.add_node(node_edge.0 as u8);
    }
    // add the edges
    for (idx, node_edge) in node_edge_list.iter().enumerate() {
        for edge in &node_edge.1 {
            graph.add_edge(NodeIndex::new(idx), NodeIndex::new(edge.0), edge.1 as i32);
        }
    }
    //graph.add_node(weight);
    // show graph 
    //println!("{}", Dot::new(&graph.map(|_, n| (*n) as char, |_, e| *e)));
    graph
}

fn concancate_files () {
    let mut output = File::create("result/ml_file").unwrap();
    let inputs = vec!["data/0_mldata.txt", "data/1_mldata.txt", "data/2_mldata.txt", "data/3_mldata.txt", "data/4_mldata.txt", "data/5_mldata.txt", "data/6_mldata.txt", "data/7_mldata.txt", "data/9_mldata.txt", "data/10_mldata.txt", "data/12_mldata.txt"];
    for i in inputs {
        let mut input = File::open(i).unwrap();
        io::copy(&mut input, &mut output).unwrap();
    }
    println!("done");
}

pub fn get_quality_score_count_confident_error () {
    let mut quality_score_count: Vec<usize> = vec![0; 94];
    // get the error locations
    let error_locations = get_error_bases_from_himut_vcf (); //chromosone, location, ref allele, alt allele
    let confident_locations = get_confident_locations_from_file ();
    for confident_location in &confident_locations {
        for error_location in &error_locations {
            let error_chromosone_result = error_location.0[3..].parse::<usize>();
            let error_chromosone_num;
            match error_chromosone_result {
                Ok(x) => {error_chromosone_num = x;},
                Err(_) => {continue;},
            }
            let confident_chromosone_num = confident_location.0[3..].parse::<usize>().unwrap();
            if error_chromosone_num > confident_chromosone_num {
                break;
            }
            else if error_chromosone_num < confident_chromosone_num {
                continue;
            }
            else {
                // condition that the error and confident have same chr number
                if confident_location.2 < error_location.1 {
                    // error location has passed confident location
                    break;
                }
                if confident_location.1 > error_location.1 {
                    // error location is behind the location
                    continue;
                }
                if confident_location.1 < error_location.1 && confident_location.2 > error_location.1 {
                    // within range
                    let path = READ_BAM_PATH;
                    let mut bam_reader = BamIndexedReader::from_path(path).unwrap();
                    let temp_quality_scores = get_quality_scores_and_base_at_location (error_location.1, 1, &error_location.0, &mut bam_reader);
                    for quality_score in temp_quality_scores {
                        if quality_score.1 as char == error_location.3 {
                            quality_score_count[quality_score.0 as usize] += 1;
                            println!("chr{}: {} < {} < {}, {} -> {} ({})", confident_chromosone_num, confident_location.1, error_location.1, confident_location.2, error_location.2, quality_score.1 as char, quality_score.0);
                        }
                    }
                }
            }
        }
    }
    println!("{:?}", quality_score_count);
}

pub fn get_quality_score_count_confident (thread_id: usize) {
    let mut quality_score_count: Vec<usize> = vec![0; 94];
    // read the bed file find the confident regions
    let confident_locations = get_confident_locations_from_file ();
    let chr_to_process = format!("chr{}", thread_id + 1);
    for confident_location in confident_locations {
        if confident_location.0 != chr_to_process {
            continue;
        }
        println!("{:?}", confident_location);
        let path = READ_BAM_PATH;
        let mut bam_reader = BamIndexedReader::from_path(path).unwrap();
        let temp_quality_scores = get_quality_scores_and_base_at_location (confident_location.1, confident_location.2 - confident_location.1, &confident_location.0, &mut bam_reader);
        for quality_score in temp_quality_scores {
            quality_score_count[quality_score.0 as usize] += 1;
        }
    }
    println!("{:?}", quality_score_count);
    let write_file = format!("result/{}_confidentquality.txt", thread_id);
    let write_string = format!("{:?}", quality_score_count);
    write_string_to_file(&write_file, &write_string);
}

fn get_confident_locations_from_file () -> Vec<(String, usize, usize)> {
    let mut location_vec: Vec<(String, usize, usize)> = vec![];
    let file_path = CONFIDENT_PATH;
    let f = File::open(&file_path).unwrap();
    let mut reader = BufReader::new(f);
    let mut buffer = String::new();
    loop {
        buffer.clear();
        match reader.read_line(&mut buffer) {
            Ok(_) => {
                let mut split_text_iter = (buffer.split("\t")).into_iter();
                let chromosone ; 
                match split_text_iter.next() {
                    Some(x) => {chromosone = x.to_string();},
                    None => {break;},
                }
                let start_string;
                match split_text_iter.next() {
                    Some(x) => {start_string = x.to_string();},
                    None => {break;},
                }
                let mut end_string;
                match split_text_iter.next() {
                    Some(x) => {end_string = x.to_string();},
                    None => {break;},
                };
                end_string.pop();
                let end = end_string.parse::<usize>().unwrap();
                let start = start_string.parse::<usize>().unwrap();
                //println!("{} {} {}",chromosone, start, end);
                location_vec.push((chromosone, start, end));
            },
            Err(_) => {break;},
        };
    }
    location_vec
}

pub fn get_data_for_ml (start: usize, end: usize, thread_id: usize) {
    let chromosone = format!("{}{}", String::from("chr"), 1);
    let mut position_base = start;
    let mut error_index = 0;
    let error_locations = get_error_bases_from_himut_vcf (); //chromosone, location, ref allele, alt allele
    // get the error index of required chromosone
    loop {
        if error_locations[error_index].0.eq(&chromosone) && error_locations[error_index].1 > start {
            break;
        }
        error_index += 1;
    }
    'bigloop: loop {
        if position_base % 1000 == 0 {
            println!("Thread ID: {} Position {}", thread_id, position_base);
        }
        let seq_name_qual_and_errorpos_vec = get_corrosponding_seq_name_location_quality_from_bam(position_base, &chromosone.to_string(), &'X');
        for seq_name_qual_and_errorpos in &seq_name_qual_and_errorpos_vec {
            // get the three base context
            let mut fai_reader = faidx::Reader::from_path(&"/data1/GiaB_benchmark/GRCh38.fa").unwrap();
            let threebase_context = read_fai_file(position_base - 1, &chromosone, &mut fai_reader);

            let mut error = false;
            // error is here
            if error_locations[error_index].1 == position_base {
                if seq_name_qual_and_errorpos.0.as_bytes()[seq_name_qual_and_errorpos.3] == error_locations[error_index].3 as u8 {
                    error = true;
                }
            }
            let base = seq_name_qual_and_errorpos.0.as_bytes()[seq_name_qual_and_errorpos.3] as char;
            let quality = seq_name_qual_and_errorpos.2;

            let parallel_stuff;
            // check if the file is already available
            if check_file_availability(&seq_name_qual_and_errorpos.1, INTERMEDIATE_PATH) {
                let available_file_path = format!("{}/{}", INTERMEDIATE_PATH, seq_name_qual_and_errorpos.1);
                let temp_quality_score = get_quality_scores_from_file(&available_file_path, seq_name_qual_and_errorpos.3);
                parallel_stuff = temp_quality_score.2;
            }
            else {
                continue;
            }
            // write data
            let write_string = format!("{} {} {} {} {}", error, threebase_context, base, quality, parallel_stuff);
            let write_file = format!("result/{}_mldata.txt", thread_id);
            write_string_to_file(&write_file, &write_string);
        }
        if error_locations[error_index].1 == position_base {
            error_index += 1;
        }
        
        // get the ccs at the location
        position_base += 1;
        if position_base > end {
            break 'bigloop;
        }
    }

}

pub fn get_quality_score_count_topology_cut_errors (start: usize, end: usize, thread_id: usize) {
    let mut quality_score_count: Vec<usize> = vec![0; 94];
    let chromosone = format!("{}{}", String::from("chr"), 21);
    let error_locations = get_error_bases_from_himut_vcf (); //chromosone, location, ref allele, alt allele
    // go through the error locations
    let mut index = 0;
    for error_location in error_locations {
        if (error_location.1 < start) || !(error_location.0.eq(&chromosone)) {
            continue;
        }
        if (error_location.1 > end) && (error_location.0.eq(&chromosone)) {
            break;
        }
        index += 1;
        let seq_name_qual_and_errorpos_vec = get_corrosponding_seq_name_location_quality_from_bam(error_location.1, &chromosone.to_string(), &'X');
        for seq_name_qual_and_errorpos in &seq_name_qual_and_errorpos_vec {
            // check if the file is already available
            if check_file_availability(&seq_name_qual_and_errorpos.1, INTERMEDIATE_PATH) {
                let available_file_path = format!("{}/{}", INTERMEDIATE_PATH, seq_name_qual_and_errorpos.1);
                let temp_quality_score = get_quality_scores_from_file(&available_file_path, seq_name_qual_and_errorpos.3);
                if error_location.3 == temp_quality_score.0 as char {
                    quality_score_count[temp_quality_score.1 as usize] += 1;
                }
            }
        }
    }
    println!("count {}", index);
    let write_file = format!("result/{}_errorcount.txt", thread_id);
    let write_string = format!("{:#?}", quality_score_count);
    write_string_to_file(&write_file, &write_string);
}

pub fn get_quality_score_count_topology_cut (start: usize, end: usize, thread_id: usize) {
    let mut quality_score_count: Vec<usize> = vec![0; 200];
    let chromosone = format!("{}{}", String::from("chr"), 21);
    let mut position_base = start; 
    'bigloop: loop {
        if position_base % 1000 == 0 {
            println!("Thread ID: {} Position {}", thread_id, position_base);
        }
        // get the required info from sam
        let seq_name_qual_and_errorpos_vec = get_corrosponding_seq_name_location_quality_from_bam(position_base, &chromosone.to_string(), &'X');
        for seq_name_qual_and_errorpos in &seq_name_qual_and_errorpos_vec {
            // check if the file is already available
            if check_file_availability(&seq_name_qual_and_errorpos.1, INTERMEDIATE_PATH) {
                let available_file_path = format!("{}/{}", INTERMEDIATE_PATH, seq_name_qual_and_errorpos.1);
                let temp_quality_score = get_quality_scores_from_file(&available_file_path, seq_name_qual_and_errorpos.3);
                quality_score_count[temp_quality_score.1 as usize] += 1;
            }
            if position_base > end {
                break 'bigloop;
            }
        }
        position_base += 1;
    }
    let write_file = format!("result/{}_totalcount.txt", thread_id);
    let write_string = format!("{:#?}", quality_score_count);
    write_string_to_file(&write_file, &write_string);
}

pub fn get_quality_scores_from_file(file_path: &String, required_pos: usize) -> (u8, u8, String) {
    let mut temp_quality_vec: (u8, u8, String) = (0, 0, "".to_string());
    // open the file
    let f = File::open(&file_path).unwrap();
    let mut reader = BufReader::new(f);
    let mut buffer = String::new();
    let mut current_pos = 0;
    loop {
        buffer.clear();
        reader.read_line(&mut buffer).unwrap();
        if current_pos == required_pos {
            let mut split_text_iter = (buffer.split(" ")).into_iter();
            let base;
            match split_text_iter.next() {
                Some(x) => {if x.as_bytes().len() > 0 {base = x.as_bytes()[0];} else {break;}},
                None => {break;},
            }
            let temp;
            match split_text_iter.next() {
                Some(x) => {temp = x.parse::<u8>();},
                None => {break;}
            }
            let quality;
            match temp {
                Ok(x) => {quality = x;},
                Err(_) => {break;}
            }
            //println!("{} {}", base, quality);
            temp_quality_vec = (base, quality, buffer);
            break;
        }
        current_pos += 1;
    }
    temp_quality_vec
}

pub fn pipeline_process_all_ccs_file_poa (chromosone: &str, start: usize, end: usize, thread_id: usize) {
    let mut index_thread = 0;
    let mut skip_thousand = false;
    let mut skip_index = 0;
    'bigloop: for process_location in start..end {
        // skip thousand when same found
        if skip_thousand {
            skip_index += 1;
            if skip_index > 10000 {
                skip_thousand = false;
                skip_index = 0;
            }
            else {
                continue;
            }
        }
        println!("Thread {}: Chr {} Loc {}, tasks_done {}", thread_id, chromosone, process_location, index_thread);
        // get the string and the name
        let seq_name_qual_and_errorpos_vec = get_corrosponding_seq_name_location_quality_from_bam(process_location, &chromosone.to_string(), &'X');
        let mut all_skipped = true;
        for seq_name_qual_and_errorpos in &seq_name_qual_and_errorpos_vec {
            println!("Thread {}: Processing ccs file: {}", thread_id, seq_name_qual_and_errorpos.1);
            // check if the file is already available
            if check_file_availability(&seq_name_qual_and_errorpos.1, INTERMEDIATE_PATH) {
                println!("Thread {}: File Available, skipping", thread_id);
                continue;
            }
            // if not available do poa and make a file
            // find the subreads of that ccs
            let mut sub_reads = get_the_subreads_by_name_sam(&seq_name_qual_and_errorpos.1);
            // skip if no subreads, errors and stuff
            if sub_reads.len() == 0 {
                continue;
            }
            all_skipped = false;
            // filter out the long reads and rearrange the reads
            sub_reads = reverse_complement_filter_and_rearrange_subreads(&sub_reads);
            // reverse if score is too low
            sub_reads = check_the_scores_and_change_alignment(sub_reads, &seq_name_qual_and_errorpos.0);

            sub_reads.insert(0, seq_name_qual_and_errorpos.0.clone());
            // do poa with the read and subreads, get the poa and consensus
            let mut sequence_number: usize = 0;
            let mut aligner = Aligner::new(MATCH, MISMATCH, GAP_OPEN, &sub_reads[0].as_bytes().to_vec(), BAND_SIZE);
            
            for sub_read in &sub_reads {
                if sequence_number != 0 {
                    aligner.global(&sub_read.as_bytes().to_vec()).add_to_graph();
                }
                let node_num = aligner.graph().node_count();
                if node_num > MAX_NODES_IN_POA {
                    println!("NUM OF NODES {} TOO BIG, SKIPPING", node_num);
                    skip_thousand = true;
                    continue 'bigloop;
                }
                sequence_number += 1;
                println!("Thread {}: Sequence {} processed", thread_id, sequence_number);
                if sequence_number > 10 {
                    break;
                }
            }
            let (calculated_consensus, calculated_topology) = aligner.poa.consensus(); //just poa
            let calculated_graph: &Graph<u8, i32, Directed, usize> = aligner.graph();
            let quality_output = get_consensus_quality_scores(sequence_number, &calculated_consensus, &calculated_topology, calculated_graph);
            // match the calculated consensus to the original consensus and get the required indices
            let calculated_indices = get_redone_consensus_matched_positions(&sub_reads[0], &calculated_consensus);
            let mut index = 0;
            for byte in sub_reads[0].as_bytes() {
                let character = *byte as char;
                let calculated_index = calculated_indices[index];
                let write_string = format!("{} {} {} {:?}", character, quality_output.0[calculated_index] as usize, (sub_reads.len() - 1) ,quality_output.1[calculated_index]);
                let write_file = format!("{}/{}", INTERMEDIATE_PATH, &seq_name_qual_and_errorpos.1);
                write_string_to_file(&write_file, &write_string);
                index += 1;
            }
            index_thread += 1;
        }
        if all_skipped {
            skip_thousand = true;
        }
    }
}

fn check_file_availability (file_name: &str, search_path: &str) -> bool {
    let temp_path_string = format!("{}/{}", search_path, file_name);
    let path = std::path::Path::new(&temp_path_string);
    let prefix = path.parent().unwrap();
    let file_available = match read_dir(prefix) {
        Ok(_) => {true},
        Err(_) => {false},
    };
    file_available
}

fn get_redone_consensus_matched_positions (pacbio_consensus: &String, calculated_consensus: &Vec<u8>) -> Vec<usize> {
    let mut consensus_matched_indices: Vec<usize> = vec![];
    let pacbio_consensus_vec: Vec<u8> = pacbio_consensus.bytes().collect();
    let (alignment, _) = pairwise(&calculated_consensus, &pacbio_consensus_vec, MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND);
    let mut calc_index = 0;
    for op in alignment {
        match op as char {
            'm' => {
                consensus_matched_indices.push(calc_index);
                calc_index += 1;
            },
            's' => {
                calc_index += 1;
                consensus_matched_indices.push(0);
            },
            'd' => {
                consensus_matched_indices.push(0);
            },
            'i' => {
                calc_index += 1;
            },
            _ => {},
        }
    }
    consensus_matched_indices
}

pub fn pipeline_redo_poa_get_topological_quality_score (chromosone: &str, start: usize, end: usize, thread_id: usize) {
    // get the error locations
    let error_locations = get_error_bases_from_himut_vcf (); //chromosone, location, ref allele, alt allele
    // go through the error locations
    let mut index = 0;
    let time_instant = Instant::now();
    for error_location in error_locations {
        if (error_location.0 == chromosone) && (error_location.1 < start) {
            continue;
        }
        if (error_location.0 != chromosone) && (error_location.1 > end) {
            break;
        }
        println!("Thread ID: {}, Error position {}:{} ref allele: {} alt allele: {}", thread_id, error_location.0, error_location.1, error_location.2, error_location.3);
        // find the ccs which are in that error
        let seq_name_qual_and_errorpos_vec = get_corrosponding_seq_name_location_quality_from_bam(error_location.1, &error_location.0, &error_location.3);
        for seq_name_qual_and_errorpos in seq_name_qual_and_errorpos_vec {
            println!("Processing ccs file: {}", seq_name_qual_and_errorpos.1);
            // find the subreads of that ccs
            let mut sub_reads = get_the_subreads_by_name_sam(&seq_name_qual_and_errorpos.1);
            // skip if no subreads, errors and stuff
            if sub_reads.len() == 0 {
                continue;
            }
            // filter out the long reads and rearrange the reads
            sub_reads = reverse_complement_filter_and_rearrange_subreads(&sub_reads);
            // reverse the sub reads if score is low, very low probability
            //sub_reads = check_the_scores_and_change_alignment(sub_reads, &seq_name_qual_and_errorpos.0);
            
            sub_reads.insert(0, seq_name_qual_and_errorpos.0.clone());
            println!("CURRENT BAND SIZE = {}", BAND_SIZE);
            // do poa with the read and subreads, get the poa and consensus
            let mut sequence_number: usize = 0;
            let mut aligner = Aligner::new(MATCH, MISMATCH, GAP_OPEN, &sub_reads[0].as_bytes().to_vec(), BAND_SIZE);
            
            for sub_read in &sub_reads {
                if sequence_number != 0 {
                    aligner.global(&sub_read.as_bytes().to_vec()).add_to_graph();
                }
                sequence_number += 1;
                println!("Sequence {} processed", sequence_number);
            }
            let (calculated_consensus, calculated_topology) = aligner.poa.consensus(); //just poa
            let calculated_graph: &Graph<u8, i32, Directed, usize> = aligner.graph();
            // check if fixed (check the consensus location which is matched to the read error location)
            let pacbio_error_pos_node_index = seq_name_qual_and_errorpos.3;
            println!("BASE {} ", seq_name_qual_and_errorpos.0.as_bytes()[pacbio_error_pos_node_index] as char);
            let position = get_redone_consensus_error_position(&seq_name_qual_and_errorpos.0, &calculated_consensus, seq_name_qual_and_errorpos.3);
            //match calculated_topology.iter().position(|r| *r == seq_name_qual_and_errorpos.3) {
            //    Some(x) => {position = x;},
            //    None => {position = get_redone_consensus_error_position(&seq_name_qual_and_errorpos.0, &calculated_consensus, seq_name_qual_and_errorpos.3);},
            //}
            println!("pacbio position {} calculated position {},", seq_name_qual_and_errorpos.3, position);
            // calculate the quality score of the location
            let skip_nodes: Vec<usize> = calculated_topology[0 .. position + 1].to_vec();
            let target_node_parent;
            let target_node_child;
            if position == 0 {
                target_node_parent = None;
            }
            else {
                target_node_parent = Some(calculated_topology[position - 1]);
            }
            if (position + 1) >= calculated_consensus.len() {
                target_node_child = None;
            }
            else {
                target_node_child = Some(calculated_topology[position + 1]);
            }
            let (parallel_nodes, parallel_num_incoming_seq, _) = get_parallel_nodes_with_topology_cut (skip_nodes, sequence_number,  calculated_topology[position], target_node_parent, target_node_child, calculated_graph);
            let (calculated_quality_score, _, parallel_bases, _) = base_quality_score_calculation (sequence_number, parallel_nodes, parallel_num_incoming_seq, calculated_consensus[position], calculated_graph);
            let write_string = format!("Error position {}:{} ref allele: {} alt allele: {}\nPacbio base: \t{} quality: {}\nCalculated base: \t{} quality: {}\nParallel Bases: ACGT:{:?}\n\n", error_location.0, error_location.1, error_location.2, error_location.3, error_location.3, seq_name_qual_and_errorpos.2, calculated_consensus[position] as char, calculated_quality_score, parallel_bases);
            let write_file = format!("result/quality_{}.txt", thread_id);
            write_string_to_file(&write_file, &write_string);
            let write_string = format!("{}\n{}\n\n", write_string, get_zoomed_graph_section(calculated_graph, &calculated_topology[position]));
            let write_file = format!("result/graph_{}.txt", thread_id);
            write_string_to_file(&write_file, &write_string);
            index += 1;
            println!("Thread {} task done, total time:{:?} total tasks:{}", thread_id, time_instant.elapsed() ,index);
        }
    }
}

fn reverse_complement_filter_and_rearrange_subreads (original_subreads: &Vec<String>) -> Vec<String> {
    let mut seqvec: Vec<String> = vec![];
    //reverse complement every other line
    let mut index = 0;
    for seq in original_subreads {
        if index % 2 != 0 {
            let mut tempseq: Vec<char> = vec![];
            let iterator = seq.chars().rev().into_iter();
            for char in iterator{
                tempseq.push(match char {
                    'A' => 'T',
                    'C' => 'G',
                    'G' => 'C',
                    'T' => 'A',
                    _ => ' ',
                });
            }
            seqvec.push(tempseq.iter().cloned().collect::<String>());
        }
        else {
            seqvec.push((*seq.clone()).to_string());
        }
        index += 1;
    }
    //get rid of the last incomplete reading
    //seqvec.pop();
    //sort the vector by size
    seqvec.sort_by_key(|seq| seq.len());
    //drop the sequences which are > 1.8x median size
    let median_size: f32 = seqvec[(seqvec.len() / 2) - 1].len() as f32;
    if seqvec.len() > 5 {
        let mut drop_index = seqvec.len();
        for index in (seqvec.len() / 2)..(seqvec.len() - 1) {
            if seqvec[index].len() as f32 > (median_size * 1.5) {
                drop_index = index;
                break;
            }
        }
        for _ in drop_index..seqvec.len() {
            seqvec.pop();
        }
    }
    // rearrange the seq vector median first and rest according mediand size difference
    seqvec.sort_by(|a, b| ((a.len() as f32 - median_size).abs()).partial_cmp(&(b.len() as f32 - median_size).abs()).unwrap());
    seqvec
}

fn get_the_subreads_by_name_sam (full_name: &String) -> Vec<String> {
    let mut subread_vec: Vec<String> = vec![];
    let mut split_text_iter = (full_name.split("/")).into_iter();
    let file_name = split_text_iter.next().unwrap();
    let required_id = split_text_iter.next().unwrap().parse::<usize>().unwrap();
    let path = format!("{}{}{}", DATA_PATH.to_string(), file_name, ".subreads.sam".to_string());
    if file_name.eq(&"m64125_201017_124255".to_string()) {
        return subread_vec;
    }
    let file_position = read_index_file_for_sam (&file_name.to_string(), required_id);
    // file stuff init
    let mut count = 0;
    let f = File::open(&path).unwrap();
    let mut reader = BufReader::new(f);
    let mut buffer = String::new();
    reader.seek(SeekFrom::Start(file_position as u64)).expect("");
    loop {
        buffer.clear();
        reader.read_line(&mut buffer).unwrap();
        // split it to find the id
        let mut temp_split_iter = (buffer.split("/")).into_iter();
        temp_split_iter.next();
        let current_id;
        match temp_split_iter.next().unwrap().parse::<usize>() {
            Ok(x) => {current_id = x;},
            Err(_) => {break;},
        }
        if current_id != required_id {
            break;
        }
        else {
            // write code to extract the sequence and add to subread_vec
            let mut data_split_iter = (buffer.split("\t")).into_iter();
            println!("{}", data_split_iter.next().unwrap());
            for _ in 0..8 {data_split_iter.next();}
            subread_vec.push(data_split_iter.next().unwrap().to_string());
            count += 1;
        }
    }
    println!("count = {}", count);
    drop(reader);
    drop(buffer);
    subread_vec
}

pub fn read_index_file_for_sam (file_name: &String, read_name: usize) -> usize {
    // get the file location from index file if no index found make one
    println!("Reading index file {}{}", file_name, ".cui".to_string());
    let index_path = format!("result/{}.cui", file_name);
    let f;
    f = match File::open(&index_path) {
        Ok(x) => {x},
        Err(_) => {make_index_file_for_sam(&file_name.to_string())},
    };
    let mut reader = BufReader::new(f);
    let mut buffer = String::new();
    
    // get the file length
    reader.seek(SeekFrom::End(0)).expect("");
    let end_file_pos = reader.stream_position().unwrap();
    
    //go back to start
    reader.seek(SeekFrom::Start(0)).expect("");

    // go through and find index
    // jump to the middle of the file **binary searching**
    let mut section_start_file_pos = 0;
    let mut section_end_file_pos = end_file_pos;
    let mut current_file_pos = (section_start_file_pos + section_end_file_pos) / 2;
    let mut required_position = 0;
    loop {
        reader.seek(SeekFrom::Start(current_file_pos)).expect("");
        // get rid of the half line
        buffer.clear();
        reader.read_line(&mut buffer).unwrap();
        // the required line
        buffer.clear();
        reader.read_line(&mut buffer).unwrap();
        // split it to find the id
        let mut temp_split_iter = (buffer.split("\t")).into_iter();
        
        let current_id;
        match temp_split_iter.next().unwrap().parse::<usize>() {
            Ok(x) => {current_id = x;},
            Err(_) => {break;},
        }
        let required_position_string = temp_split_iter.next().unwrap().replace("\n", "");
        required_position = required_position_string.parse::<usize>().unwrap();
        //println!("curr: {}", current_id);
        // jumping 
        if read_name == current_id {
            break;
        }
        else if read_name > current_id {
            section_start_file_pos = current_file_pos;
        }
        else {
            section_end_file_pos = current_file_pos;
        }
        current_file_pos = (section_start_file_pos + section_end_file_pos) / 2;
    }
    required_position
}

pub fn make_index_file_for_sam (file_name: &String) -> File {
    println!("Making index file for {}{}", file_name, ".subreads.sam".to_string());
    let path = format!("{}{}{}", DATA_PATH.to_string(), file_name, ".subreads.sam".to_string());
    let write_path = format!("result/{}.cui", file_name);
    // get the file name and load it
    // file stuff init
    let f = File::open(&path).unwrap();
    let mut reader = BufReader::new(f);
    let mut buffer = String::new();

    // get the file length
    reader.seek(SeekFrom::End(0)).expect("");
    let end_file_pos = reader.stream_position().unwrap();
    
    //go back to start
    reader.seek(SeekFrom::Start(0)).expect("");
    
    // go through the file saving the sequence indices
    let mut write_string: String = "".to_string();
    let mut current_position;
    let mut current_ccs: usize;
    let mut prev_ccs: usize = 0;
    let mut index = 0;
    let mut break_count = 0;
    loop {
        // get the next line if available
        buffer.clear();
        match reader.read_line(&mut buffer) {
            Ok(_) => {if buffer.len() == 0 {break;}},
            Err(_) => {break;},
        }
        // split it to find the id
        let mut temp_split_iter = (buffer.split("/")).into_iter();
        temp_split_iter.next();
        match temp_split_iter.next() {
            Some(x) => {current_ccs = match x.parse::<usize>() { Ok(val) => {val}, Err(_) => {0}};},
            None => {continue;},
        }
        // add to the pos_read_name if different
        if current_ccs != prev_ccs {
            break_count = 0;
            prev_ccs = current_ccs;
            current_position = reader.stream_position().unwrap();
            write_string = format!("{}\n{}\t{}", write_string, current_ccs, current_position);
            index += 1;
            // display progress and write the current data
            if index % 100000 == 0 {
                println!("Progress {}%", (current_position * 100) / end_file_pos);
                write_string_to_file(&write_path, &write_string);
                write_string = "".to_string();
            }
        }
        else {
            break_count += 1;
            if break_count > 10000 {
                break;
            }
        }
    }
    // write the rest
    write_string_to_file(&write_path, &write_string);
    File::open(&write_path).unwrap()
}

fn get_the_subreads_by_name_bam (error_chr: &String, error_pos: usize, full_name: &String) -> Vec<String> {
    let subread_vec: Vec<String> = vec![];
    let mut split_text_iter = (full_name.split("/")).into_iter();
    let file_name = split_text_iter.next().unwrap();
    let required_id = split_text_iter.next().unwrap().parse::<i64>().unwrap();
    let path = format!("{}{}{}", DATA_PATH.to_string(), file_name, ".subreads.mapped.bam".to_string());
    if file_name.eq(&"m64125_201017_124255".to_string()) {
        return subread_vec;
    }
    let mut bam_reader = BamIndexedReader::from_path(path).unwrap();
    bam_reader.fetch((error_chr, error_pos as i64, error_pos as i64 + 1)).unwrap();
    let mut index = 0;
    for read in bam_reader.records() {
        let readunwrapped = read.unwrap();
        // get the required name
        let read_name = String::from_utf8(readunwrapped.qname().to_vec()).expect("");
        let mut split_read_name_iter = (read_name.split("/")).into_iter();
        split_read_name_iter.next().unwrap();
        let css_id = split_read_name_iter.next().unwrap().parse::<i64>().unwrap();
        if required_id != css_id {
            continue;
        }

        println!("readname {}", read_name);
        if readunwrapped.seq_len() < 5 {
            continue;
        }
        index += 1;
    }
    println!("new count = {}", index);
    
    subread_vec
}

fn get_redone_consensus_error_position (pacbio_consensus: &String, calculated_consensus: &Vec<u8>, pacbio_error_position: usize) -> usize {
    let mut consensus_match_invalid_indices: Vec<usize> = vec![];
    let pacbio_consensus_vec: Vec<u8> = pacbio_consensus.bytes().collect();
    let mut aligned_pacbio_scores_vec: Vec<usize> = vec![];
    let mut aligned_pacbio_bases:Vec<u8> = vec![];
    let (alignment, _) = pairwise(&calculated_consensus, &pacbio_consensus_vec, MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND);
    let mut pacbio_index = 0;
    let mut calc_index = 0;
    let mut calc_error_position: usize = 0;
    for op in alignment {
        match op as char {
            'm' => {
                pacbio_index += 1;
                calc_index += 1;
            },
            's' => {
                aligned_pacbio_bases.push(pacbio_consensus_vec[pacbio_index]);
                consensus_match_invalid_indices.push(calc_index);
                pacbio_index += 1;
                calc_index += 1;
            },
            'd' => {
                pacbio_index += 1;
            },
            'i' => {
                consensus_match_invalid_indices.push(calc_index);
                aligned_pacbio_bases.push(126);
                aligned_pacbio_scores_vec.push(33);
                calc_index += 1;
            },
            _ => {},
        }
        if pacbio_error_position == pacbio_index {
            calc_error_position = calc_index
        }
    }
    println!("Pacbio error position {} corrosponds to calculated error position {}", pacbio_error_position, calc_error_position);
    calc_error_position
}

fn check_the_scores_and_change_alignment (seqvec: Vec<String>, pacbio_consensus: &String) -> Vec<String> {
    let mut invert: bool = false;
    let mut seqvec2: Vec<String> = vec![];
    // check the scores for 3 sequences
    let mut index = 0;
    for seq in &seqvec {
        let (_, score) = pairwise(&pacbio_consensus.as_bytes().to_vec(), &seq.as_bytes().to_vec(), MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND);
        println!("score: {}", score);
        if score < 2000 {
            invert = true;
            break;
        }
        else if index > 1 {
            break;
        }
        index += 1;
    }
    if invert {
        println!("Scores are too low, inverting sequences.");
        //reverse complement every line
        for seq in &seqvec {
            let mut tempseq: Vec<char> = vec![];
            let iterator = seq.chars().rev().into_iter();
            for char in iterator{
                tempseq.push(match char {
                    'A' => 'T',
                    'C' => 'G',
                    'G' => 'C',
                    'T' => 'A',
                    _ => ' ',
                });
            }
            seqvec2.push(tempseq.iter().cloned().collect::<String>());
        }
    }
    else {
        seqvec2 = seqvec;
    }
    seqvec2
}

fn get_corrosponding_seq_name_location_quality_from_bam (error_pos: usize, error_chr: &String, base_change: &char) -> Vec<(String, String, u8, usize)> {
    let mut seq_name_qual_and_errorpos: Vec<(String, String, u8, usize)> = vec![];
    let path = &READ_BAM_PATH;
    let mut bam_reader = BamIndexedReader::from_path(path).unwrap();
    bam_reader.fetch((error_chr, error_pos as i64, error_pos as i64 + 1)).unwrap();
    'read_loop: for read in bam_reader.records() {
        let readunwrapped = read.unwrap();
        // get the data
        let mut read_index = 0;
        let read_name = String::from_utf8(readunwrapped.qname().to_vec()).expect("");
        let read_vec = readunwrapped.seq().as_bytes().to_vec();
        let read_string = String::from_utf8(readunwrapped.seq().as_bytes().to_vec()).expect("");
        if readunwrapped.seq_len() < 5 {
            continue;
        }
        // get the location from the cigar processing
        let mut temp_character_vec: Vec<char> = vec![];
        // get the read start position
        let read_start_pos = readunwrapped.pos() as usize;
        let mut current_ref_pos = read_start_pos;
        let mut current_read_pos = 0;
        // decode the cigar string
        for character in readunwrapped.cigar().to_string().as_bytes() {
            match *character as char {
                'M' => {         
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    if (current_ref_pos + temp_int >= error_pos)
                        && (current_ref_pos <= error_pos + 1) {
                        (read_index, _) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, error_pos, 1);
                        if &(read_vec[read_index] as char) == base_change && ('X' != *base_change) {
                            break;
                        }
                        else if 'X' == *base_change {
                            break;
                        }
                        else {
                            continue 'read_loop;
                        }
                    }
                    current_ref_pos += temp_int;
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'H' => {
                    temp_character_vec = vec![];
                },
                'S' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'I' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'N' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    if (current_ref_pos + temp_int >= error_pos)
                        && (current_ref_pos <= error_pos + 1) {
                        //(_, _) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, error_pos, 1);
                        continue 'read_loop;
                    }
                    current_ref_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'D' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    if (current_ref_pos + temp_int >= error_pos)
                        && (current_ref_pos <= error_pos + 1) {
                        //let (_, _) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, error_pos, 1);
                        continue 'read_loop;
                    }
                    current_ref_pos += temp_int;
                    temp_character_vec = vec![];
                },
                _ => {
                    temp_character_vec.push(*character as char);
                },
            }
        }
        seq_name_qual_and_errorpos.push((read_string.clone(), read_name.clone(), readunwrapped.qual()[read_index], read_index));
    }
    drop(bam_reader);
    seq_name_qual_and_errorpos
}

pub fn write_string_to_file (file_name: &String, input_string: &String) {
    let path = std::path::Path::new(&file_name);
    let prefix = path.parent().unwrap();
    create_dir_all(prefix).unwrap();
    let mut file = OpenOptions::new().create(true).append(true).open(file_name).unwrap();
    writeln!(file, "{}", input_string).expect("result file cannot be written");
}

pub fn get_zoomed_graph_section (normal_graph: &Graph<u8, i32, Directed, usize>, focus_node: &usize)-> String {
    let mut graph_section= "".to_string();
    let normal_dot = format!("{:?}", Dot::new(&normal_graph.map(|_, n| (*n) as char, |_, e| *e)));
    let displaying_nodes: Vec<usize> = find_neighbouring_indices (NUM_OF_ITER_FOR_ZOOMED_GRAPHS, *focus_node, normal_graph);
    let mut graph_section_nodes: String = "".to_string();
    let mut graph_section_edges: String = "".to_string();
    //find the position in the dot file and add to the graph section
    for node in &displaying_nodes {
        //get the nodes from the dot file
        match normal_dot.find(&format!(" {} [", node)) {
            Some(start) => {
                let mut char_seq = vec![];
                let mut end = start;
                loop {
                    char_seq.push(normal_dot.chars().nth(end).unwrap());
                    if normal_dot.chars().nth(end).unwrap() == '\n' {
                        break;
                    }
                    end += 1;
                }
                //add from start to end to the graph_section string
                graph_section_nodes = char_seq.iter().collect::<String>();
            },
            None => {}
        }
        graph_section = format!("{}{}", graph_section, graph_section_nodes).to_string();
        graph_section_nodes = "".to_string();
    }
    for node in &displaying_nodes {
        //get the edges from the dot file
        let edge_entries: Vec<usize> = normal_dot.match_indices(&format!(" {} ->", node)).map(|(i, _)|i).collect();
        for edge in edge_entries {
            let mut char_seq = vec![];
                let mut end = edge;
                loop {
                    char_seq.push(normal_dot.chars().nth(end).unwrap());
                    if normal_dot.chars().nth(end).unwrap() == '\n' {
                        break;
                    }
                    end += 1;
                }
                //add from start to end to the graph_section string
                graph_section_edges = format!("{}{}", graph_section_edges, char_seq.iter().collect::<String>()).to_string();
        }
        graph_section = format!("{}{}", graph_section, graph_section_edges).to_string();
        graph_section_edges = "".to_string();
    }
    //modifying the section graph with the highlight
    graph_section = modify_dot_graph_with_highlight (graph_section, focus_node);
    //make it a dot graph
    graph_section = format!("digraph {{\n{} }}", graph_section);
    graph_section
}


fn find_neighbouring_indices (num_of_iterations: usize, focus_node: usize, graph: &Graph<u8, i32, Directed, usize> ) -> Vec<usize> {
    let mut indices: Vec<usize> = vec![];
    if num_of_iterations <= 0 {
        return indices;
    }
    let mut immediate_neighbours = graph.neighbors_undirected(NodeIndex::new(focus_node));
    while let Some(neighbour_node) = immediate_neighbours.next() {
        if !indices.contains(&neighbour_node.index()){
            indices.push(neighbour_node.index());
        }
        let obtained_indices = find_neighbouring_indices(num_of_iterations - 1, neighbour_node.index(), graph);
        for obtained_index in obtained_indices {
            if !indices.contains(&obtained_index){
                indices.push(obtained_index);
            }
        }
    }
    indices 
}

fn modify_dot_graph_with_highlight (mut dot: String, focus_node: &usize) -> String {
    match dot.find(&format!(" {} [", focus_node)) {
        Some(mut x) => {
            while dot.chars().nth(x).unwrap() != '[' {
                x += 1;
            }
            dot.replace_range(x + 12..x + 12,"Target node");
        },
        None => {}
    };
    dot
}

pub fn get_error_bases_from_himut_vcf () -> Vec<(String, usize, char, char)> {
    let mut error_locus_vec: Vec<(String, usize, char, char)> = vec![]; //chromosone, position, ref, alt
    let path = &"data/somatic.vcf";
    let mut bcf = Reader::from_path(path).expect("Error opening file.");
    // iterate through each row of the vcf body.
    for (_, record_result) in bcf.records().enumerate() {
        let record = record_result.expect("Fail to read record");
        let mut allele_vec: Vec<char> = vec![];
        for allele in record.alleles() {
            for c in allele {
                allele_vec.push(char::from(*c));
            }
        }
        let temp_str = record.desc();
        let mut split_text_iter = (temp_str.split(":")).into_iter();
        let chromosone = split_text_iter.next().unwrap();
        error_locus_vec.push((chromosone.to_string(), record.pos() as usize, allele_vec[0], allele_vec[1]));
    }
    println!("number of errors = {}", error_locus_vec.len());
    error_locus_vec
}

fn get_error_quality_score_count (error_locus_vec: Vec<(String, usize, char, char)>) {
    let mut quality_score_count: Vec<usize> = vec![0; 94];
    // read the merged mapped sorted bam file
    let path = &"data/merged.bam";
    let mut bam_reader = BamIndexedReader::from_path(path).unwrap();
    let mut index = 0;
    // go through the errors and update the count
    for error_locus in error_locus_vec {
        let position_base = error_locus.1;
        let temp_quality_scores_and_bases = get_quality_scores_and_base_at_location (position_base, 1, &error_locus.0, &mut bam_reader);
        for (quality_score, base) in temp_quality_scores_and_bases {
            if error_locus.3 == base as char {
                quality_score_count[quality_score as usize] += 1;
            }
        }
        if index % 1000 == 0 {
            println!("currently processing error {}", index);
        }
        index += 1;
    }
    println!("{:#?}", quality_score_count);
}

fn get_quality_scores_and_base_at_location (required_pos: usize, required_len: usize, chromosone: &String, reader: &mut BamIndexedReader) -> Vec<(u8, u8)> {
    // reused function, this one will only work for a single position
    let mut quals: Vec<(u8, u8)> = vec![]; //quality score, base
    
    match reader.fetch((chromosone, required_pos as i64, required_pos as i64 + required_len as i64)) {
        Ok(_) => {},
        Err(_) => {return vec![]},
    }
    for read in reader.records() {
        let readunwrapped = read.unwrap();
        if readunwrapped.seq_len() < required_len {
            continue;
        }
        let mut temp_character_vec: Vec<char> = vec![];
        // get the read start position
        let read_start_pos = readunwrapped.pos() as usize;
        let mut current_ref_pos = read_start_pos;
        let mut current_read_pos = 0;
        // decode the cigar string
        for character in readunwrapped.cigar().to_string().as_bytes() {
            match *character as char {
                'M' => {         
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    if (current_ref_pos + temp_int >= required_pos)
                        && (current_ref_pos <= required_pos + required_len) {
                        let (p1, p2) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, required_pos, required_len);
                        for index in p1..p2 {
                            quals.push((readunwrapped.qual()[index], readunwrapped.seq().as_bytes()[p1]));
                        }
                    }
                    current_ref_pos += temp_int;
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'H' => {
                    temp_character_vec = vec![];
                },
                'S' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'I' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'N' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    current_ref_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'D' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    current_ref_pos += temp_int;
                    temp_character_vec = vec![];
                },
                _ => {
                    temp_character_vec.push(*character as char);
                },
            }
        }
    }
    quals
}

fn pipeline_3base_context () {
    // make a vector with 256 entries for each correct and incorrect count eg entry index AAA -> A = 0 & TTT -> T = 255
    let mut count_vector: Vec<usize> = vec![0; 256];
    // go from chr1 to chr21
    for index in 1..2 {
        let chromosone = format!("{}{}", String::from("chr"), index.to_string());
        println!("Reading {}", chromosone);
        let mut position_base = 1090000;
        // go from 1mil to 24mil bases in small lengths
        loop {
            if position_base % 100000 == 0 {
                println!("Position {}", position_base);
            }
            //get the reads and reference
            let path = &"data/merged.bam";
            let mut bam_reader = BamIndexedReader::from_path(path).unwrap();
            let path = &"data/GRCh38.fa";
            let mut fai_reader = faidx::Reader::from_path(path).unwrap();
            let reads = read_bam_file(position_base, &chromosone, &mut bam_reader);
            let reference = read_fai_file(position_base, &chromosone, &mut fai_reader);
            // process
            find_error_in_three_base_context(&reference, &reads, &mut count_vector);
            position_base += THREE_BASE_CONTEXT_READ_LENGTH;
            if position_base > 5000000 {
                break;
            }
        }
    }
    // print the results
    print_3base_context_results(&count_vector);
}

fn read_bam_file (required_start_pos: usize, chromosone: &String, reader: &mut BamIndexedReader) -> Vec<String> {
    let mut read_vec: Vec<String> = vec![];
    reader.fetch((chromosone, required_start_pos as i64, required_start_pos as i64 + THREE_BASE_CONTEXT_READ_LENGTH as i64 - 1)).unwrap();
    for read in reader.records() {
        let readunwrapped = read.unwrap();
        if readunwrapped.seq_len() < THREE_BASE_CONTEXT_READ_LENGTH {
            continue;
        }
        let mut temp_character_vec: Vec<char> = vec![];
        let mut temp_read_vec: Vec<u8> = vec![];
        // get the read start position
        let read_start_pos = readunwrapped.pos() as usize;
        let mut current_ref_pos = read_start_pos;
        let mut current_read_pos = 0;
        // decode the cigar string
        for character in readunwrapped.cigar().to_string().as_bytes() {
            match *character as char {
                'M' => {         
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    //println!("M RUNNIN at pos {} for {}", current_ref_pos, temp_int); 
                    if (current_ref_pos + temp_int >= required_start_pos)
                        && (current_ref_pos <= required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                        let (p1, p2) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, required_start_pos, THREE_BASE_CONTEXT_READ_LENGTH);
                        temp_read_vec = [temp_read_vec, readunwrapped.seq().as_bytes()[p1..p2].to_vec()].concat();
                    }
                    current_ref_pos += temp_int;
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'H' => {
                    temp_character_vec = vec![];
                },
                'S' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'I' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'N' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    //println!("N RUNNIN at pos {} for {}", current_ref_pos, temp_int); 
                    if (current_ref_pos + temp_int >= required_start_pos)
                        && (current_ref_pos <= required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                        let (p1, p2) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, required_start_pos, THREE_BASE_CONTEXT_READ_LENGTH);
                        for _ in p1..p2 {
                            temp_read_vec.push('X' as u8);
                        }
                    }
                    current_ref_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'D' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    //println!("D RUNNIN at pos {} for {}", current_ref_pos, temp_int); 
                    if (current_ref_pos + temp_int >= required_start_pos)
                        && (current_ref_pos <= required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                        let (p1, p2) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, required_start_pos, THREE_BASE_CONTEXT_READ_LENGTH);
                        for _ in p1..p2 {
                            temp_read_vec.push('X' as u8);
                        }
                    }
                    current_ref_pos += temp_int;
                    temp_character_vec = vec![];
                },
                _ => {
                    temp_character_vec.push(*character as char);
                },
            }
        }
        // take care of cases when string is not long enough to cover the whole ref required
        read_vec.push(String::from_utf8(temp_read_vec).expect("Found invalid UTF-8"));
    }
    read_vec
}

fn read_fai_file (start_pos: usize, chromosone: &String, reader: &mut faidx::Reader) -> String {
    //println!("The ref sequenece");
    //println!("{}", fasta_reader.fetch_seq_string(chromosone, start_pos, start_pos + THREE_BASE_CONTEXT_READ_LENGTH).unwrap());
    reader.fetch_seq_string(chromosone, start_pos, start_pos + THREE_BASE_CONTEXT_READ_LENGTH).unwrap()
}

pub fn get_required_start_end_positions_from_read (section_length: usize, current_ref_pos: usize, current_read_pos: usize, required_pos: usize, required_len: usize) -> (usize, usize) {
    let p1;
    let p2;
    // case when both end partial match is required
    if (current_ref_pos < required_pos) && (current_ref_pos + section_length > required_pos + required_len) {
        p1 = required_pos - current_ref_pos + current_read_pos;
        p2 = current_read_pos + required_pos + required_len - current_ref_pos;
    }
    // case when start partial matches are required
    else if (current_ref_pos < required_pos) && (current_ref_pos + section_length >= required_pos) {
        p1 = required_pos - current_ref_pos + current_read_pos;
        p2 = current_read_pos + section_length;
    }
    // case when end partial matches are required
    else if (current_ref_pos <= required_pos + required_len) && (current_ref_pos + section_length > required_pos + required_len) {
        p1 = current_read_pos;
        p2 = current_read_pos + required_pos + required_len - current_ref_pos;
    }
    // normal case
    else {
        p1 = current_read_pos;
        p2 = current_read_pos + section_length;
    }
    (p1, p2)
}

fn homopolymer_test () {
    // get random generated data
    let seqvec = get_random_sequences_from_generator(RANDOM_SEQUENCE_LENGTH, NUMBER_OF_RANDOM_SEQUENCES, SEED);
    // run the rustbio poa
    println!("Processing seq 1");
    let mut aligner = Aligner::new(MATCH, MISMATCH, GAP_OPEN, &seqvec[0].as_bytes().to_vec(), 100);
    let mut index = 0;
    for seq in &seqvec {
        if index != 0 {
            println!("Processing seq {}", index + 1);
            aligner.global(&seq.as_bytes().to_vec()).add_to_graph();
        }
        index += 1;
    }
    let (consensus, _) = aligner.poa.consensus();
    let consensus_score = get_consensus_score(&seqvec, &consensus, MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND);
    println!("rustbio score = {}", consensus_score);
    for base in &consensus {
        print!("{}", *base as char);
    }
    println!("");
    // run the homopolymer poa
    // convert the sequence to homopolymer
    let mut homopolymervec: Vec<Vec<HomopolymerCell>> = vec![];
    for seq in &seqvec {
        homopolymervec.push(convert_sequence_to_homopolymer(seq));
    }
    println!("Processing seq 1");
    let mut aligner = Poa::initialize(&homopolymervec[0], MATCH, MISMATCH, GAP_OPEN);
    let mut index = 0;
    for seq in &homopolymervec {
        if index != 0 {
            println!("Processing seq {}", index + 1);
            aligner.add_to_poa(&seq);
        }
        index += 1;
    }
    let (consensus, _consensus_node_indices) = aligner.consensus();
    let consensus_score = get_consensus_score(&seqvec, &consensus, MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND);
    //println!("{}", format!("{:?}", Dot::new(&aligner.poa_graph.map(|_, n| (*n) as char, |_, e| *e))));
    println!("homopolymer score = {}", consensus_score);
    for base in &consensus {
        print!("{}", *base as char);
    }
    println!("");

}

pub fn find_error_in_three_base_context (reference: &String, reads: &Vec<String>, count_vector: &mut Vec<usize>) {
    // ignore the edges
    let mut index = 1;
    loop {
        // check if index is at the end of the ref
        if index + 1 >= reference.len() {
            break;
        } 
        // get the three base around index in reference
        let ref_3base: Vec<u8> = reference.as_bytes()[index - 1..index + 2].to_vec();
        // print the ref 3 bases
        //print!("ref base: ");
        for _base in &ref_3base {
            //print!("{}", *base as char);
        }
        //println!("");
        // get the bits of ref_3base
        let ref_6bit = get_6bit_from_3base (&ref_3base);
        // go through the reads and count the errors and correct ones
        for read in reads {
            if !(index + 1 >= read.len()) {
                let read_base: u8 = read.as_bytes()[index];
                //println!("read base: {}", read_base as char);
                // get the bits of read_3base
                let error_base_bit = get_2bit_from_base (&read_base);
                // concancate the two together
                let count_vector_index = (ref_6bit * 4 + error_base_bit) as usize;
                // add the entry to the count_vector
                count_vector[count_vector_index] += 1;
            }
        }
        index += 1;
    }
    //println!("{:?}", count_vector);
    //print_3base_context_results (&count_vector);
}

pub fn print_3base_context_results (count_vector: &Vec<usize>) {
    let mut stats_for_3bases: Vec<(usize, usize)> = vec![(0, 0); 64];
    let mut checker = 3;
    let mut temp_correct_count = 0;
    let mut temp_wrong_count = 0;
    for index in 0..256 {
        let (current_base, mutation, correct) = get_3base_mutation_result_from_8bit (index as u8);
        for base in &current_base {
            print!("{}", *base as char);
        }
        print!(" -> {}, count: {}", mutation as char, count_vector[index]);
        println!("");
        if correct {
            temp_correct_count += count_vector[index];
        }
        else {
            temp_wrong_count += count_vector[index];
        }
        if index >= checker {
            stats_for_3bases[(index - 3) / 4] = (temp_correct_count, temp_wrong_count);
            temp_correct_count = 0;
            temp_wrong_count = 0;
            checker += 4;
        }
    }
    println!("STATS");
    for index in 0..64 {
        let percentage = stats_for_3bases[index].1 as f64 / (stats_for_3bases[index].1  + stats_for_3bases[index].0) as f64;
        for base in get_3base_from_6bit(index as u8) {
            print!("{}", base as char);
        }
        print!(" Correct: {:>8} Error: {:>8} Error rate: {:>8}", stats_for_3bases[index].0, stats_for_3bases[index].1, percentage);
        println!("");
    }
}
pub fn get_3base_from_6bit (_6bit: u8) -> Vec<u8> {
    let mut _3base = vec![];
    let mut mask = 3;
    // make the _3base vector
    let mut index = 0;
    let mut multiplier = 1;
    loop {
        _3base.push(get_base_from_2bit((_6bit & mask) / multiplier));
        index += 1;
        if index == 3 {
            break;
        }
        mask *= 4;
        multiplier *= 4;
    }
    _3base.reverse();
    _3base
}

pub fn get_6bit_from_3base (_3base: &Vec<u8>) -> u8 {
    assert!(_3base.len() == 3);
    let mut _3base_bit = 0;
    let mut multiplier: u8 = 16;
    for base in _3base {
        match *base as char {
            'A' => {_3base_bit += 0 * multiplier},
            'C' => {_3base_bit += 1 * multiplier},
            'G' => {_3base_bit += 2 * multiplier},
            'T' => {_3base_bit += 3 * multiplier},
            _ => {},
        }
        multiplier /= 4;
    } 
    return _3base_bit;
}

pub fn get_3base_mutation_result_from_8bit (_8bit: u8) -> (Vec<u8>, u8, bool) {
    // mask the last two bits {AND operator with 3 = 00000011} and get mutation
    let mutation = get_base_from_2bit (_8bit & 3);
    // get the last two 
    let _3base_bit = _8bit / 4;
    // make the _3base vector
    let _3base = get_3base_from_6bit (_3base_bit);
    let correct = if _3base[1] == mutation {
        true
    }
    else {
        false
    };
    (_3base, mutation, correct)
}

pub fn get_base_from_2bit (_2bit: u8) -> u8 {
    match _2bit {
        0 => {'A' as u8},
        1 => {'C' as u8},
        2 => {'G' as u8},
        3 => {'T' as u8},
        _ => {'X' as u8},
    }
}

pub fn get_2bit_from_base (base: &u8) -> u8 {
    match *base as char {
        'A' => {0},
        'C' => {1},
        'G' => {2},
        'T' => {3},
        _ => {0},
    }
}

pub fn convert_sequence_to_homopolymer (sequence: &String) -> Vec<HomopolymerCell> {
    let mut homopolymer_vec: Vec<HomopolymerCell> = vec![];
    let mut prev_base: u8 = 0;
    let mut frequency: usize = 1;
    for base in sequence.as_bytes().to_vec() {
        if prev_base == base {
            frequency += 1;
        }
        else if prev_base != 0 {
            homopolymer_vec.push(HomopolymerCell::new(prev_base, frequency));
            frequency = 1;
        }
        prev_base = base;
    }
    homopolymer_vec.push(HomopolymerCell::new(prev_base, frequency));
    homopolymer_vec
}

pub fn get_consensus_score(seqvec : &Vec<String>, consensus: &Vec<u8>, match_score: i32, mismatch_score: i32, gap_open_score: i32, gap_extend_score: i32) -> isize {
    let mut consensus_score = 0;
    for seq in seqvec {
        let (_, score) = pairwise(&consensus, &seq.as_bytes().to_vec(), match_score, mismatch_score, gap_open_score, gap_extend_score);
        consensus_score += score;
    }
    consensus_score
}

#[derive(Clone)]
pub struct HomopolymerCell {
    pub base: u8,
    pub frequency: usize,
}

impl HomopolymerCell {
    pub fn new(base: u8, frequency: usize) -> Self {
        HomopolymerCell {
            base: base,
            frequency: frequency,
        }
    }
}