mod alignment;
mod generator;

use alignment::pairwise::pairwise;
use generator::simple::get_random_sequences_from_generator;
use alignment::poa::Poa;
use petgraph::dot::Dot;

const SEED: u64 = 2;
const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = -2;
const MATCH: i32 = 2;
const MISMATCH: i32 = -3;
const RANDOM_SEQUENCE_LENGTH: usize = 200;
const NUMBER_OF_RANDOM_SEQUENCES: usize = 5;

fn main() {
    let seqvec = get_random_sequences_from_generator(RANDOM_SEQUENCE_LENGTH, NUMBER_OF_RANDOM_SEQUENCES, SEED);
    let (_, score) = pairwise(&seqvec[0].as_bytes().to_vec(), &seqvec[1].as_bytes().to_vec(), MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND);
    println!("Pairwise score = {}", score);
    println!("Processing seq 1");
    let mut aligner = Poa::initialize(&seqvec[0].as_bytes().to_vec(), MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND);
    let mut index = 0;
    for seq in &seqvec {
        if index != 0 {
            println!("Processing seq {}", index + 1);
            aligner.add_to_poa(&seq.as_bytes().to_vec());
        }
        index += 1;
    }
    let (consensus, _consensus_node_indices) = aligner.consensus();
    let consensus_score = get_consensus_score(&seqvec, &consensus);
    println!("{}", format!("{:?}", Dot::new(&aligner.poa_graph.map(|_, n| (*n) as char, |_, e| *e))));
    println!("consensus score == {}", consensus_score);
    for base in &consensus {
        print!("{}", *base as char);
    }
    println!("");
}

fn get_consensus_score(seqvec : &Vec<String>, consensus: &Vec<u8>) -> isize {
    let mut consensus_score = 0;
    for seq in seqvec {
        let (_, score) = pairwise(&consensus, &seq.as_bytes().to_vec(), MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND);
        consensus_score += score;
    }
    consensus_score
}