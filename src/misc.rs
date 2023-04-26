use crate::alignment::pairwise::pairwise;
#[derive(Clone)]
pub struct HomopolymerCell {
    pub base: u8,
    pub frequency: usize,
}

impl HomopolymerCell {
    fn new(base: u8, frequency: usize) -> Self {
        HomopolymerCell {
            base: base,
            frequency: frequency,
        }
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