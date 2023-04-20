mod alignment;
mod generator;

use alignment::pairwise::pairwise;
use generator::simple::get_random_sequences_from_generator;
use alignment::poa::Poa;

const SEED: u64 = 0;
const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = -1;
const MATCH: i32 = 2;
const MISMATCH: i32 = -2;
const RANDOM_SEQUENCE_LENGTH: usize = 10;
const NUMBER_OF_RANDOM_SEQUENCES: usize = 2;

fn main() {
    let test = get_random_sequences_from_generator(RANDOM_SEQUENCE_LENGTH, NUMBER_OF_RANDOM_SEQUENCES, SEED);
    let (_, score) = pairwise(&test[0].as_bytes().to_vec(), &test[1].as_bytes().to_vec(), MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND);
    println!("score = {}", score);
    let test1 = Poa::initialize(&test[0].as_bytes().to_vec(), MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND);
    Poa::get_alignment(test1, &test[1].as_bytes().to_vec());
    
}
