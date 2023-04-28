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
use misc::HomopolymerCell;
use rust_htslib::bcf::{Reader, Read};
use std::convert::TryFrom;

const SEED: u64 = 2;
const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = 0;
const MATCH: i32 = 2;
const MISMATCH: i32 = -2;
const RANDOM_SEQUENCE_LENGTH: usize = 100;
const NUMBER_OF_RANDOM_SEQUENCES: usize = 5;

fn main() {
    homopolymer_test();
    //read_vcf_file();
}

fn read_vcf_file () {
    let path = &"data/output.vcf.gz";
    let mut bcf = Reader::from_path(path).expect("Error opening file.");
    // iterate through each row of the vcf body.
    for (_i, record_result) in bcf.records().enumerate() {
        let record = record_result.expect("Fail to read record");
        let mut s = String::new();
        for allele in record.alleles() {
            for c in allele {
                s.push(char::from(*c))
            }
            s.push(' ')
        }
        // 0-based position and the list of alleles
        println!("Locus: {}, Alleles: {}", record.pos(), s);
        // number of sample in the vcf
        let sample_count = usize::try_from(record.sample_count()).unwrap();
        // Counting ref, alt and missing alleles for each sample
        let mut n_ref = vec![0; sample_count];
        let mut n_alt = vec![0; sample_count];
        let mut n_missing = vec![0; sample_count];
        let gts = record.genotypes().expect("Error reading genotypes");
        for sample_index in 0..sample_count {
            // for each sample
            for gta in gts.get(sample_index).iter() {
                // for each allele
                match gta.index() {
                    Some(0) => n_ref[sample_index] += 1,  // reference allele
                    Some(_) => n_alt[sample_index] += 1,  // alt allele
                    None => n_missing[sample_index] += 1, // missing allele
                }
            }
        }
    }
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