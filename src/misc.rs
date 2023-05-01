use crate::alignment::pairwise::pairwise;
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

pub fn find_error_in_three_base_context (reference: &String, reads: &Vec<String>) {
    assert!(reference.len() >= 3);
    for read in reads {
        assert!(read.len() >= 3);
    }
    // make a vector with 4096 entries for each correct and incorrect count eg entry index AAA -> AAA = 0 & TTT -> TTT = 4095
    let mut count_vector: Vec<(usize, usize)> = vec![(0, 0); 4096];
    // ignore the edges
    let mut index = 1;
    loop {
        // check if index is at the end of the ref
        if index + 1 >= reference.len() {
            break;
        } 
        // get the three base around index in reference
        let mut ref_3base: Vec<u8> = reference.as_bytes()[index - 1..index + 2].to_vec();
        // get the bits of ref_3base

        // go through the reads and count the errors and correct ones
        for read in reads {
            if !(index + 1 >= read.len()) {
                let read_3base: Vec<u8> = read.as_bytes()[index - 1..index + 2].to_vec();
                // check if the middle base is same or not
                let correct = check_if_middle_base_same (&ref_3base, &read_3base);
                // get the bits of read_3base

                // concancate the two together

                let count_vector_index = 0;
                // add the entry to the count_vector
                if correct {
                    count_vector[count_vector_index].1 += 1;
                }
                else {
                    count_vector[count_vector_index].0 += 1;
                }
            }
        }
        index += 1;
    }
    


}

pub fn check_if_middle_base_same (ref_3base: &Vec<u8>, read_3base: &Vec<u8>) -> bool {
    assert!(ref_3base.len() == 3);
    assert!(read_3base.len() == 3);
    if ref_3base[1] == read_3base[1] {
        true
    }
    else {
        false
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