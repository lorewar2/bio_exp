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
use std::convert::TryFrom;

const SEED: u64 = 2;
const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = 0;
const MATCH: i32 = 2;
const MISMATCH: i32 = -2;
const RANDOM_SEQUENCE_LENGTH: usize = 1000;
const NUMBER_OF_RANDOM_SEQUENCES: usize = 5;

fn main() {
    //homopolymer_test();
    //let seqvec = get_random_sequences_from_generator(RANDOM_SEQUENCE_LENGTH, NUMBER_OF_RANDOM_SEQUENCES, SEED);
    //find_error_in_three_base_context (&seqvec[0], &seqvec);
    //read_fai_file();
    read_bam_file();
}

fn read_bam_file () {
    let path = &"/Users/wmw0016/Documents/mount1/data1/hifi_consensus/try2/merged.bam";
    let mut bam = IndexedReader::from_path(path).unwrap();
    bam.fetch(("chr1", 10005, 10500)).unwrap(); // coordinates 10000..20000 on reference named "chrX"
    let mut index = 0;
    for read in bam.records() {
        //println!("read name: {:?}", read.unwrap().qname());
        let readunwrapped = read.unwrap();
        println!("read {} {:?}", readunwrapped.pos(), readunwrapped.seq().as_bytes());
        for character in readunwrapped.qname() {
            print!("{}", *character as char);
        }
        break;
        println!("");
        index += 1;
    }
    println!("{}", index);
}

fn read_fai_file () {
    let path = &"data/GRCh38.fa";
    let mut fasta_reader = faidx::Reader::from_path(path).expect("Error opening file.");
    println!("{}", fasta_reader.fetch_seq_string("chr1", 10000, 10500).unwrap());
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