

#![allow(dead_code)]
//#[global_allocator]
//static ALLOC: dhat::Alloc = dhat::Alloc;
mod alignment;
mod generator;
mod misc;
mod quality;
//use pprof;
//use crate::misc::pipeline_redo_poa_get_topological_quality_score;
//use crate::misc::pipeline_process_all_ccs_file_poa;
//use std::thread;
use crate::misc::get_quality_score_count_confident_error;
//use crate::misc::get_quality_score_count_topology_cut;
//use crate::misc::get_data_for_ml;
//use crate::alignment::poabandedsmarter::Aligner;
//use crate::generator::simple::get_random_sequences_from_generator;
//use crate::alignment::pairwise::pairwise;
const SEED: u64 = 2;
const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = 0;
const MATCH: i32 = 2;
const MISMATCH: i32 = -2;
const NTHREADS: usize = 8;

fn main() {
    get_quality_score_count_confident_error ();
    // Make a vector to hold thfe children which are spawneddd.
    //let _profiler = dhat::Profiler::new_heap();
    /* 
    let mut children = vec![];
    let chromosone = "chr1";
    let total_start = 10_000_000;
    let total_end = 250_000_000;
    let one_thread_allocation = (total_end - total_start) / NTHREADS;
    for i in 0..NTHREADS {
        // calculate my start and end locations
        // Spin up another thread
        children.push(thread::spawn(move || {
            // calculate my start and end locations
            let start = total_start + one_thread_allocation * i;
            let end = start + one_thread_allocation;
            println!("Thread number {} started, {} from {} to {}..", chromosone, i, start, end);
            pipeline_process_all_ccs_file_poa (chromosone, start, end, i);
            //get_quality_score_count_topology_cut(start, end, i);
            //get_data_for_ml (start, end, i);
            //pipeline_redo_poa_get_topological_quality_score(chromosone, start, end, i);
        }));
    }

    for child in children {
        // Wait for the thread to finish. Returns a result.
        let _ = child.join();
    }
    //
     */
    /*let sequences = get_random_sequences_from_generator(2000, 2, 6);
    let mut sequence_number = 0;
    let mut aligner = Aligner::new(MATCH, MISMATCH, GAP_OPEN, &sequences[0].as_bytes().to_vec(), 20);
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


