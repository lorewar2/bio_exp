#![allow(dead_code)]
mod alignment;
mod generator;
mod misc;
mod quality;
//use pprof;
use crate::generator::simple::get_random_sequences_from_generator;
use crate::alignment::poabandedsmarter::Aligner;
use crate::alignment::pairwise::pairwise;
use crate::misc::pipeline_redo_poa_get_topological_quality_score;
const SEED: u64 = 2;
const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = 0;
const MATCH: i32 = 2;
const MISMATCH: i32 = -2;

fn main() {
    pipeline_redo_poa_get_topological_quality_score();
    /* 
    let sequences = get_random_sequences_from_generator(10000, 10, 6);
    let mut sequence_number = 0;
    let mut aligner = Aligner::new(MATCH, MISMATCH, GAP_OPEN, &sequences[0].as_bytes().to_vec(), 100);
    for sequence in &sequences {
        if sequence_number != 0 {
            aligner.global(&sequence.as_bytes().to_vec()).add_to_graph();
        }
        sequence_number += 1;
        println!("Sequence {} processed", sequence_number);
    }
    let (calculated_consensus, _) = aligner.poa.consensus(); //just poa
    let mut total_quality = 0;
    for sequence in &sequences {
        let (_, temp_score) = pairwise(&calculated_consensus, &sequence.as_bytes().to_vec(), MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND);
        total_quality += temp_score;
    }
    println!("quality {}", total_quality);
    */
}


