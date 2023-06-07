use crate::alignment::pairwise::pairwise;
use crate::generator::simple::get_random_sequences_from_generator;
use crate::alignment::poarustbio::Aligner;
use crate::alignment::poahomopolymer::Poa;
use petgraph::{Graph, Directed, graph::NodeIndex};
use petgraph::dot::Dot;
use rust_htslib::bam::{Read as BamRead, IndexedReader as BamIndexedReader};
use rust_htslib::bcf::{Reader, Read as BcfRead};
use rust_htslib::faidx;
use std::{fs::OpenOptions, io::{prelude::*}, path::Path};

const SEED: u64 = 2;
const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = 0;
const MATCH: i32 = 2;
const MISMATCH: i32 = -2;
const RANDOM_SEQUENCE_LENGTH: usize = 1000;
const NUMBER_OF_RANDOM_SEQUENCES: usize = 5;
const THREE_BASE_CONTEXT_READ_LENGTH: usize = 1000;
const NUM_OF_ITER_FOR_ZOOMED_GRAPHS: usize = 4;

pub fn write_string_to_file (file_name: impl AsRef<Path>, input_string: String) {
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

fn pipeline_quality_score_error_graph () {
    // get the quality scores in the ccs
    get_quality_score_count ();
    // get the errors from himut vcf (no somatic mutations in this file)
    let error_locations = get_error_bases_from_himut_vcf();
    // get the quality scores of error positions
    get_error_quality_score_count (error_locations);
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

fn get_quality_score_count () {
    let skip_length = 3000;
    let continue_threshold = 10000;
    let mut continue_count = 0;
    let mut quality_score_count: Vec<usize> = vec![0; 94];
    // read the merged mapped sorted bam file
    let path = &"data/merged.bam";
    let mut bam_reader = BamIndexedReader::from_path(path).unwrap();
    // go from chr1 to chr21
    for index in 1..22 {
        let chromosone = format!("{}{}", String::from("chr"), index.to_string());
        println!("Reading {}", chromosone);
        let mut position_base = 5000000;
        let mut prev_93_count = usize::MAX; 
        // go from 1mil to 240mil bases in small lengths (skip length)
        loop {
            if position_base % 1000000 == 0 {
                println!("Position {}", position_base);
            }
            // iterate through by counting the quality scores.
            let temp_quality_scores = get_quality_scores_and_base_at_location (position_base, skip_length, &chromosone, &mut bam_reader);
            for quality_score in temp_quality_scores {
                quality_score_count[quality_score.0 as usize] += 1;
            }
            //println!("{:?}", quality_score_count);
            if quality_score_count[93] == prev_93_count {
                continue_count += 1;
                if continue_count >= continue_threshold {
                    break;
                }
            }
            else {
                continue_count = 0;
            }
            prev_93_count = quality_score_count[93];
            position_base += skip_length;
        }
    }
    println!("{:#?}", quality_score_count);
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
    let mut aligner = Aligner::new(MATCH, MISMATCH, GAP_OPEN, &seqvec[0].as_bytes().to_vec());
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