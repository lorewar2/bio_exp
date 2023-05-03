#![allow(dead_code)]
mod alignment;
mod generator;
mod misc;
mod quality;

use generator::simple::get_random_sequences_from_generator;
use alignment::poarustbio::Aligner;
use alignment::poahomopolymer::Poa;
//use petgraph::dot::Dot;
use misc::get_consensus_score;
use misc::convert_sequence_to_homopolymer;
use misc::find_error_in_three_base_context;
use misc::HomopolymerCell;
use rust_htslib::faidx;
use rust_htslib::bam::{Read, IndexedReader};

const SEED: u64 = 2;
const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = 0;
const MATCH: i32 = 2;
const MISMATCH: i32 = -2;
const RANDOM_SEQUENCE_LENGTH: usize = 1000;
const NUMBER_OF_RANDOM_SEQUENCES: usize = 5;
const THREE_BASE_CONTEXT_READ_LENGTH: usize = 100;

fn main() {
    //homopolymer_test();
    //let seqvec = get_random_sequences_from_generator(RANDOM_SEQUENCE_LENGTH, NUMBER_OF_RANDOM_SEQUENCES, SEED);
    //find_error_in_three_base_context (&seqvec[0], &seqvec);
    //read_fai_file();
    //
    pipeline_3base_context();
}

fn pipeline_3base_context () {
    // go from chr1 to chr21
    let mut position_base = 1000000;
    // go from 1mil to 24mil bases in small lengths
    // do 1 iteration to test
    // get the reference genome of length

    //get the reads
    read_bam_file(position_base);
    read_fai_file(position_base);

}

fn read_bam_file (required_start_pos: usize) {
    let mut read_vec: Vec<String> = vec![];
    let path = &"data/merged.bam";
    let mut bam = IndexedReader::from_path(path).unwrap();
    bam.fetch(("chr1", required_start_pos as i64, required_start_pos as i64 + THREE_BASE_CONTEXT_READ_LENGTH as i64 - 1)).unwrap();
    let mut index = 0;
    for read in bam.records() {
        index += 1;
        if index < 5 {
            continue;
        }
        let readunwrapped = read.unwrap();
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
                    println!("M RUNNIN at pos {} for {}", current_ref_pos, temp_int); 
                    if (current_ref_pos + temp_int >= required_start_pos)
                        && (current_ref_pos <= required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                        let p1;
                        let p2;
                        // case when both end partial match is required
                        if (current_ref_pos < required_start_pos) && (current_ref_pos + temp_int > required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                            p1 = required_start_pos - current_ref_pos + current_read_pos;
                            p2 = current_read_pos + required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH - current_ref_pos;;
                        }
                        // case when start partial matches are required
                        else if (current_ref_pos < required_start_pos) && (current_ref_pos + temp_int >= required_start_pos) {
                            p1 = required_start_pos - current_ref_pos + current_read_pos;
                            p2 = current_read_pos + temp_int;
                        }
                        // case when end partial matches are required
                        else if (current_ref_pos <= required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) && (current_ref_pos + temp_int > required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                            
                            p1 = current_read_pos;
                            p2 = current_read_pos + required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH - current_ref_pos;
                        }
                        // normal case
                        else {
                            p1 = current_read_pos;
                            p2 = current_read_pos + temp_int;
                        }
                        temp_read_vec = [temp_read_vec, readunwrapped.seq().as_bytes()[p1..p2].to_vec()].concat();
                    }
                    current_ref_pos += temp_int;
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'H' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    current_read_pos += temp_int;
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
                    println!("N RUNNIN at pos {} for {}", current_ref_pos, temp_int); 
                    if (current_ref_pos + temp_int >= required_start_pos)
                        && (current_ref_pos <= required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                        let p1;
                        let p2;
                        // case when both end partial match is required
                        if (current_ref_pos < required_start_pos) && (current_ref_pos + temp_int > required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                            p1 = required_start_pos - current_ref_pos + current_read_pos;
                            p2 = current_read_pos + required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH - current_ref_pos;;
                        }
                        // case when start partial matches are required
                        else if (current_ref_pos < required_start_pos) && (current_ref_pos + temp_int >= required_start_pos) {
                            p1 = required_start_pos - current_ref_pos + current_read_pos;
                            p2 = current_read_pos + temp_int;
                        }
                        // case when end partial matches are required
                        else if (current_ref_pos <= required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) && (current_ref_pos + temp_int > required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                            p1 = current_read_pos;
                            p2 = current_read_pos + required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH - current_ref_pos;
                        }
                        // normal case
                        else {
                            p1 = current_read_pos;
                            p2 = current_read_pos + temp_int;
                        }
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
                    println!("D RUNNIN at pos {} for {}", current_ref_pos, temp_int); 
                    if (current_ref_pos + temp_int >= required_start_pos)
                        && (current_ref_pos <= required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                        let p1;
                        let p2;
                        // case when both end partial match is required
                        if (current_ref_pos < required_start_pos) && (current_ref_pos + temp_int > required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                            p1 = required_start_pos - current_ref_pos + current_read_pos;
                            p2 = current_read_pos + required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH - current_ref_pos;;
                        }
                        // case when start partial matches are required
                        else if (current_ref_pos < required_start_pos) && (current_ref_pos + temp_int >= required_start_pos) {
                            p1 = required_start_pos - current_ref_pos + current_read_pos;
                            p2 = current_read_pos + temp_int;
                        }
                        // case when end partial matches are required
                        else if (current_ref_pos <= required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) && (current_ref_pos + temp_int > required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                            
                            p1 = current_read_pos;
                            p2 = current_read_pos + required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH - current_ref_pos;
                        }
                        // normal case
                        else {
                            p1 = current_read_pos;
                            p2 = current_read_pos + temp_int;
                        }
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
    for read in read_vec {
        println!("{}",read);
        println!("read ln {}", read.len());
    }
}

fn read_fai_file (start_pos: usize) {
    let path = &"data/GRCh38.fa";
    let mut fasta_reader = faidx::Reader::from_path(path).expect("Error opening file.");
    println!("The ref sequenece");
    println!("{}", fasta_reader.fetch_seq_string("chr1", start_pos, start_pos + THREE_BASE_CONTEXT_READ_LENGTH).unwrap());
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