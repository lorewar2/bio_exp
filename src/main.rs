

#![allow(dead_code)]
mod alignment;
mod generator;
mod misc;
mod quality;
use std::thread;
use crate::misc::get_quality_score_count_confident_error;

const SEED: u64 = 2;
const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = 0;
const MATCH: i32 = 2;
const MISMATCH: i32 = -2;
const NTHREADS: usize = 1;

fn main() {
    //get_quality_score_count_confident_error();
    get_quality_score_count_confident_error();
    // make a vector to hold the children which are spawned.
    //new_poa_tester();
    let mut children = vec![];
    let chromosone = "chr21";
    let total_start = 15_000_000;
    let total_end = 40_000_000;
    let one_thread_allocation = (total_end - total_start) / NTHREADS;
    for i in 0..NTHREADS {
        // calculate my start and end locations
        // spin up another thread
        children.push(thread::spawn(move || {
            // calculate my start and end locations
            let start = total_start + one_thread_allocation * i;
            let end = start + one_thread_allocation;
            println!("Thread number {} started, {} from {} to {}..", chromosone, i, start, end);
            //get_data_for_ml (chromosone, start, end, i);
        }));
    }
    for child in children {
        // wait for the thread to finish. Returns a result.
        let _ = child.join();
    }
}


