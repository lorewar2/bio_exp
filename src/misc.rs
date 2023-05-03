use crate::alignment::pairwise::pairwise;

pub fn find_error_in_three_base_context (reference: &String, reads: &Vec<String> ) {
    assert!(reference.len() >= 3);
    for read in reads {
        assert!(read.len() >= 3);
    }
    // make a vector with 256 entries for each correct and incorrect count eg entry index AAA -> A = 0 & TTT -> T = 255
    let mut count_vector: Vec<usize> = vec![0; 256];
    // ignore the edges
    let mut index = 1;
    loop {
        // check if index is at the end of the ref
        if index + 1 >= reference.len() {
            break;
        } 
        // get the three base around index in reference
        let ref_3base: Vec<u8> = reference.as_bytes()[index - 1..index + 2].to_vec();
        // print the ref 3 bases
        print!("ref base: ");
        for base in &ref_3base {
            print!("{}", *base as char);
        }
        println!("");
        // get the bits of ref_3base
        let ref_6bit = get_6bit_from_3base (&ref_3base);
        // go through the reads and count the errors and correct ones
        for read in reads {
            if !(index + 1 >= read.len()) {
                let read_base: u8 = read.as_bytes()[index];
                println!("read base: {}", read_base as char);
                // get the bits of read_3base
                let error_base_bit = get_2bit_from_base (&read_base);
                // concancate the two together
                let count_vector_index = (ref_6bit * 4 + error_base_bit) as usize;
                // add the entry to the count_vector
                count_vector[count_vector_index] += 1;
            }
        }
        index += 1;
    }
    println!("{:?}", count_vector);
    print_3base_context_results (&count_vector);
}

pub fn print_3base_context_results (count_vector: &Vec<usize>) {
    let mut stats_for_3bases: Vec<(usize, usize)> = vec![(0, 0); 64];
    let mut checker = 3;
    let mut temp_correct_count = 0;
    let mut temp_wrong_count = 0;
    for index in 0..256 {
        let (current_base, mutation, correct) = get_3base_mutation_result_from_8bit (index as u8);
        println!("current base: {:?}, mutation: {}, count: {}", current_base, mutation, count_vector[index]);
        if correct {
            temp_correct_count += count_vector[index];
        }
        else {
            temp_wrong_count += count_vector[index];
        }
        if index >= checker {
            stats_for_3bases[(index - 3) / 4] = (temp_correct_count, temp_wrong_count);
            temp_correct_count = 0;
            temp_wrong_count = 0;
            checker += 4;
        }
    }
    for index in 0..64 {
        let percentage = stats_for_3bases[index].1 as f64 / (stats_for_3bases[index].1  + stats_for_3bases[index].0) as f64;
        println!("current base: {:?}, correct {} incorrect {} percentage {}", get_3base_from_6bit(index as u8), stats_for_3bases[index].0, stats_for_3bases[index].1, percentage);
    }
}
pub fn get_3base_from_6bit (_6bit: u8) -> Vec<u8> {
    let mut _3base = vec![];
    let mut mask = 3;
    // make the _3base vector
    let mut index = 0;
    let mut multiplier = 1;
    loop {
        _3base.push(get_base_from_2bit((_6bit & mask) / multiplier));
        index += 1;
        if index == 3 {
            break;
        }
        mask *= 4;
        multiplier *= 4;
    }
    _3base.reverse();
    _3base
}

pub fn get_6bit_from_3base (_3base: &Vec<u8>) -> u8 {
    assert!(_3base.len() == 3);
    let mut _3base_bit = 0;
    let mut multiplier: u8 = 16;
    for base in _3base {
        match *base as char {
            'A' => {_3base_bit += 0 * multiplier},
            'C' => {_3base_bit += 1 * multiplier},
            'G' => {_3base_bit += 2 * multiplier},
            'T' => {_3base_bit += 3 * multiplier},
            _ => {},
        }
        multiplier /= 4;
    } 
    return _3base_bit;
}

pub fn get_3base_mutation_result_from_8bit (_8bit: u8) -> (Vec<u8>, u8, bool) {
    // mask the last two bits {AND operator with 3 = 00000011} and get mutation
    let mutation = get_base_from_2bit (_8bit & 3);
    // get the last two 
    let _3base_bit = _8bit / 4;
    // make the _3base vector
    let _3base = get_3base_from_6bit (_3base_bit);
    let correct = if _3base[1] == mutation {
        true
    }
    else {
        false
    };
    (_3base, mutation, correct)
}

pub fn get_base_from_2bit (_2bit: u8) -> u8 {
    match _2bit {
        0 => {'A' as u8},
        1 => {'C' as u8},
        2 => {'G' as u8},
        3 => {'T' as u8},
        _ => {'X' as u8},
    }
}

pub fn get_2bit_from_base (base: &u8) -> u8 {
    match *base as char {
        'A' => {0},
        'C' => {1},
        'G' => {2},
        'T' => {3},
        _ => {0},
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