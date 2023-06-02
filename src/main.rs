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
use misc::print_3base_context_results;
//use quality::topology_cut::get_consensus_quality_scores;
use misc::HomopolymerCell;
use rust_htslib::faidx;
use rust_htslib::bam::{Record, Read as BamRead, IndexedReader as BamIndexedReader, Reader as BamReader};
use rust_htslib::bcf::{Reader, Read as BcfRead};

const SEED: u64 = 2;
const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = 0;
const MATCH: i32 = 2;
const MISMATCH: i32 = -2;
const RANDOM_SEQUENCE_LENGTH: usize = 1000;
const NUMBER_OF_RANDOM_SEQUENCES: usize = 5;
const THREE_BASE_CONTEXT_READ_LENGTH: usize = 1000;
pub const MAX_USIZE: usize = 858_993_459;

fn main() {
    //homopolymer_test();
    //let seqvec = get_random_sequences_from_generator(RANDOM_SEQUENCE_LENGTH, NUMBER_OF_RANDOM_SEQUENCES, SEED);
    //pipeline_quality_score();
    //pipeline_3base_context();
    //get_quality_score_count();
    //pipeline_quality_score_error_graph ();
    pipeline_redo_poa_get_topological_quality_score();
}

fn pipeline_redo_poa_get_topological_quality_score () {
    // get the error locations
    let error_locations = get_error_bases_from_himut_vcf ();
    // go through the error locations
    let mut index = 0;
    for error_location in error_locations {
        if index < 10000{
            index += 1;
            continue;
        }
        println!("error chromosone {} position {}", error_location.0, error_location.1);
        // find the ccs which are in that error
        let seq_name_and_errorpos_vec = get_corrosponding_seq_name_location_from_bam(error_location.1, &error_location.0);
        for seq_name_and_errorpos in seq_name_and_errorpos_vec {
            println!("Processing ccs file {}", seq_name_and_errorpos.1);
            // find the subreads of that ccs
            get_the_subreads_by_name(&seq_name_and_errorpos.1);
            // do poa with the subreads

            // check if fixed
            // calculate the quality score
            break;
        }
        break;
    }
}

fn get_the_subreads_by_name (full_name: &String) -> Vec<String> {
    let mut subread_vec: Vec<String> = vec![];
    let mut split_text_iter = (full_name.split("/")).into_iter();
    //let path = format!("{}{}{}", "/Users/wmw0016/Documents/mount5/data1/hifi_consensus/try2/".to_string(), split_text_iter.next().unwrap(), ".subreads.bam".to_string());
    let path = format!("{}", "data/merged.bam");
    split_text_iter.next().unwrap();
    let ccs_name = split_text_iter.next().unwrap();

    //search for the subreads with ccs name in the file
    let mut bam = BamReader::from_path(path).unwrap();
    let mut record = Record::new();
    //bam.seek(10).unwrap();
    
    // variables for finding start position of the records and breaking right after records are read
    let average_record_len: i64 = 2500000000;
    let mut index = 0;
    let mut read_skip = false;
    let mut read_done = false;
    let mut read_started = false;
    let mut read_first_record_found = false;
    let mut jumped_behind_required_section = false;
    println!("{}", ccs_name);
    let required_id_pos = ccs_name.parse::<i64>().unwrap();
    
    while let Some(r) = bam.read(&mut record) {
        index += 1;
        
        match r {
            Ok(_) => {

            },
            Err(_) => {
                println!("current offset {}", bam.tell());
                continue;
            }
        }
        let subread_name = String::from_utf8(record.qname().to_vec()).expect("");
        println!("subreadname {}", subread_name);
        println!("subreadlen {}", record.seq_len());
        let mut sub_split_text_iter = (subread_name.split("/")).into_iter();
        let parent_ccs = sub_split_text_iter.next().unwrap();
        let current_id_pos = parent_ccs.parse::<i64>().unwrap();
        println!("current offset {} id {} ", bam.tell(), current_id_pos);
        if required_id_pos != current_id_pos {
            if !read_started {
                // if behind the required pos jump forward
                if required_id_pos > current_id_pos {
                    // estimate how many record there are between
                    let jump_value = bam.tell() + (1000) * average_record_len;
                    println!("jump value {}", jump_value);
                    // jump
                    bam.seek(jump_value).unwrap();

                }
                else {
                    read_started = true;
                }
            }
            else {
                // encountered a non match after the required records were found
                if read_first_record_found {
                    break;
                }
                if (!jumped_behind_required_section) || (required_id_pos < current_id_pos) {
                    //jump behind 10?
                    let jump_value = bam.tell() - 10 * 30 * average_record_len;
                    //jump
                    bam.seek(jump_value).unwrap();
                    jumped_behind_required_section = true;
                }
            }
        }
        else {
            read_started = true;
            if read_first_record_found {
                if record.seq_len() > 0 {
                    subread_vec.push(String::from_utf8(record.seq().as_bytes().to_vec()).expect(""));
                }
            }
            else {
                if jumped_behind_required_section {
                    read_first_record_found = true;
                    if record.seq_len() > 0 {
                        subread_vec.push(String::from_utf8(record.seq().as_bytes().to_vec()).expect(""));
                    }
                }
                else {
                    //jump behind 10?
                    let jump_value = bam.tell() - 10 * 30 * average_record_len;
                    //jump
                    bam.seek(jump_value).unwrap();
                }
            }
        }
    }
    subread_vec
}

fn get_corrosponding_seq_name_location_from_bam (error_pos: usize, error_chr: &String) -> Vec<(String, String, usize)> {
    let mut seq_name_and_errorpos: Vec<(String, String, usize)> = vec![];
    let path = &"data/merged.bam";
    let mut bam_reader = BamIndexedReader::from_path(path).unwrap();
    bam_reader.fetch((error_chr, error_pos as i64, error_pos as i64 + 1)).unwrap();
    'read_loop: for read in bam_reader.records() {
        let readunwrapped = read.unwrap();
        // get the data
        let mut read_index = 0;
        let read_name = String::from_utf8(readunwrapped.qname().to_vec()).expect("");
        let read_string = String::from_utf8(readunwrapped.seq().as_bytes().to_vec()).expect("");
        if readunwrapped.seq_len() < 5 {
            continue;
        }
        // get the location from the cigar processing
        let mut temp_character_vec: Vec<char> = vec![];
        // get the read start position
        let read_start_pos = readunwrapped.pos() as usize;
        let mut current_ref_pos = read_start_pos;
        let mut current_read_pos = 0;
        // decode the cigar string
        for character in readunwrapped.cigar().to_string().as_bytes() {
            match *character as char {
                'M' => {         
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    if (current_ref_pos + temp_int >= error_pos)
                        && (current_ref_pos <= error_pos + 1) {
                        (read_index, _) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, error_pos, 1);
                        break;
                    }
                    current_ref_pos += temp_int;
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'H' => {
                    temp_character_vec = vec![];
                },
                'S' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'I' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'N' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    if (current_ref_pos + temp_int >= error_pos)
                        && (current_ref_pos <= error_pos + 1) {
                        (_, _) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, error_pos, 1);
                        continue 'read_loop;
                    }
                    current_ref_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'D' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    if (current_ref_pos + temp_int >= error_pos)
                        && (current_ref_pos <= error_pos + 1) {
                        let (_, _) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, error_pos, 1);
                        continue 'read_loop;
                    }
                    current_ref_pos += temp_int;
                    temp_character_vec = vec![];
                },
                _ => {
                    temp_character_vec.push(*character as char);
                },
            }
        }
        seq_name_and_errorpos.push((read_string.clone(), read_name.clone(), read_index));
    }
    seq_name_and_errorpos
}

fn pipeline_quality_score_error_graph () {
    // get the quality scores in the ccs
    //get_quality_score_count ();
    // get the errors from himut vcf (no somatic mutations in this file)
    let error_locations = get_error_bases_from_himut_vcf();
    // get the quality scores of error positions
    get_error_quality_score_count (error_locations);
}

fn get_error_quality_score_count (error_locus_vec: Vec<(String, usize, char, char)>) {
    let mut quality_score_count: Vec<usize> = vec![0; 94];
    // read the merged mapped sorted bam file
    let path = &"data/merged.bam";
    let mut bam_reader = BamIndexedReader::from_path(path).unwrap();
    let mut index = 0;
    // go through the errors and update the count
    for error_locus in error_locus_vec {
        let position_base = error_locus.1;
        let temp_quality_scores_and_bases = get_quality_scores_and_base_at_location (position_base, 1, &error_locus.0, &mut bam_reader);
        for (quality_score, base) in temp_quality_scores_and_bases {
            if error_locus.3 == base as char {
                quality_score_count[quality_score as usize] += 1;
            }
        }
        if index % 1000 == 0 {
            println!("currently processing error {}", index);
        }
        index += 1;
    }
    println!("{:#?}", quality_score_count);
}

fn get_error_bases_from_himut_vcf () -> Vec<(String, usize, char, char)> {
    let mut error_locus_vec: Vec<(String, usize, char, char)> = vec![]; //chromosone, position, ref, alt
    let path = &"data/somatic.vcf";
    let mut bcf = Reader::from_path(path).expect("Error opening file.");
    // iterate through each row of the vcf body.
    for (_, record_result) in bcf.records().enumerate() {
        let record = record_result.expect("Fail to read record");
        let mut allele_vec: Vec<char> = vec![];
        for allele in record.alleles() {
            for c in allele {
                allele_vec.push(char::from(*c));
            }
        }
        let temp_str = record.desc();
        let mut split_text_iter = (temp_str.split(":")).into_iter();
        let chromosone = split_text_iter.next().unwrap();
        error_locus_vec.push((chromosone.to_string(), record.pos() as usize, allele_vec[0], allele_vec[1]));
    }
    println!("number of errors = {}", error_locus_vec.len());
    error_locus_vec
}

fn get_quality_score_count () {
    let skip_length = 3000;
    let continue_threshold = 10000;
    let mut continue_count = 0;
    let mut quality_score_count: Vec<usize> = vec![0; 94];
    // read the merged mapped sorted bam file
    let path = &"data/merged.bam";
    let mut bam_reader = BamIndexedReader::from_path(path).unwrap();
    // go from chr1 to chr21
    for index in 1..22 {
        let chromosone = format!("{}{}", String::from("chr"), index.to_string());
        println!("Reading {}", chromosone);
        let mut position_base = 5000000;
        let mut prev_93_count = MAX_USIZE; 
        // go from 1mil to 240mil bases in small lengths (skip length)
        loop {
            if position_base % 1000000 == 0 {
                println!("Position {}", position_base);
            }
            // iterate through by counting the quality scores.
            let temp_quality_scores = get_quality_scores_and_base_at_location (position_base, skip_length, &chromosone, &mut bam_reader);
            for quality_score in temp_quality_scores {
                quality_score_count[quality_score.0 as usize] += 1;
            }
            //println!("{:?}", quality_score_count);
            if quality_score_count[93] == prev_93_count {
                continue_count += 1;
                if continue_count >= continue_threshold {
                    break;
                }
            }
            else {
                continue_count = 0;
            }
            prev_93_count = quality_score_count[93];
            position_base += skip_length;
        }
    }
    println!("{:#?}", quality_score_count);
}

fn get_quality_scores_and_base_at_location (required_pos: usize, required_len: usize, chromosone: &String, reader: &mut BamIndexedReader) -> Vec<(u8, u8)> {
    // reused function, this one will only work for a single position
    let mut quals: Vec<(u8, u8)> = vec![]; //quality score, base
    
    match reader.fetch((chromosone, required_pos as i64, required_pos as i64 + required_len as i64)) {
        Ok(_) => {},
        Err(_) => {return vec![]},
    }
    for read in reader.records() {
        let readunwrapped = read.unwrap();
        if readunwrapped.seq_len() < required_len {
            continue;
        }
        let mut temp_character_vec: Vec<char> = vec![];
        // get the read start position
        let read_start_pos = readunwrapped.pos() as usize;
        let mut current_ref_pos = read_start_pos;
        let mut current_read_pos = 0;
        // decode the cigar string
        for character in readunwrapped.cigar().to_string().as_bytes() {
            match *character as char {
                'M' => {         
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    if (current_ref_pos + temp_int >= required_pos)
                        && (current_ref_pos <= required_pos + required_len) {
                        let (p1, p2) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, required_pos, required_len);
                        for index in p1..p2 {
                            quals.push((readunwrapped.qual()[index], readunwrapped.seq().as_bytes()[p1]));
                        }
                    }
                    current_ref_pos += temp_int;
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'H' => {
                    temp_character_vec = vec![];
                },
                'S' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'I' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'N' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    current_ref_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'D' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    current_ref_pos += temp_int;
                    temp_character_vec = vec![];
                },
                _ => {
                    temp_character_vec.push(*character as char);
                },
            }
        }
    }
    quals
}

fn get_required_start_end_positions_from_read (section_length: usize, current_ref_pos: usize, current_read_pos: usize, required_pos: usize, required_len: usize) -> (usize, usize) {
    let p1;
    let p2;
    // case when both end partial match is required
    if (current_ref_pos < required_pos) && (current_ref_pos + section_length > required_pos + required_len) {
        p1 = required_pos - current_ref_pos + current_read_pos;
        p2 = current_read_pos + required_pos + required_len - current_ref_pos;
    }
    // case when start partial matches are required
    else if (current_ref_pos < required_pos) && (current_ref_pos + section_length >= required_pos) {
        p1 = required_pos - current_ref_pos + current_read_pos;
        p2 = current_read_pos + section_length;
    }
    // case when end partial matches are required
    else if (current_ref_pos <= required_pos + required_len) && (current_ref_pos + section_length > required_pos + required_len) {
        p1 = current_read_pos;
        p2 = current_read_pos + required_pos + required_len - current_ref_pos;
    }
    // normal case
    else {
        p1 = current_read_pos;
        p2 = current_read_pos + section_length;
    }
    (p1, p2)
}

fn pipeline_3base_context () {
    // make a vector with 256 entries for each correct and incorrect count eg entry index AAA -> A = 0 & TTT -> T = 255
    let mut count_vector: Vec<usize> = vec![0; 256];
    // go from chr1 to chr21
    for index in 1..2 {
        let chromosone = format!("{}{}", String::from("chr"), index.to_string());
        println!("Reading {}", chromosone);
        let mut position_base = 1090000;
        // go from 1mil to 24mil bases in small lengths
        loop {
            if position_base % 100000 == 0 {
                println!("Position {}", position_base);
            }
            //get the reads and reference
            let path = &"data/merged.bam";
            let mut bam_reader = BamIndexedReader::from_path(path).unwrap();
            let path = &"data/GRCh38.fa";
            let mut fai_reader = faidx::Reader::from_path(path).unwrap();
            let reads = read_bam_file(position_base, &chromosone, &mut bam_reader);
            let reference = read_fai_file(position_base, &chromosone, &mut fai_reader);
            // process
            find_error_in_three_base_context(&reference, &reads, &mut count_vector);
            position_base += THREE_BASE_CONTEXT_READ_LENGTH;
            if position_base > 5000000 {
                break;
            }
        }
    }
    // print the results
    print_3base_context_results(&count_vector);
}



fn read_bam_file (required_start_pos: usize, chromosone: &String, reader: &mut BamIndexedReader) -> Vec<String> {
    let mut read_vec: Vec<String> = vec![];
    reader.fetch((chromosone, required_start_pos as i64, required_start_pos as i64 + THREE_BASE_CONTEXT_READ_LENGTH as i64 - 1)).unwrap();
    for read in reader.records() {
        let readunwrapped = read.unwrap();
        if readunwrapped.seq_len() < THREE_BASE_CONTEXT_READ_LENGTH {
            continue;
        }
        let mut temp_character_vec: Vec<char> = vec![];
        let mut temp_read_vec: Vec<u8> = vec![];
        // get the read start position
        let read_start_pos = readunwrapped.pos() as usize;
        let mut current_ref_pos = read_start_pos;
        let mut current_read_pos = 0;
        // decode the cigar string
        for character in readunwrapped.cigar().to_string().as_bytes() {
            match *character as char {
                'M' => {         
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    //println!("M RUNNIN at pos {} for {}", current_ref_pos, temp_int); 
                    if (current_ref_pos + temp_int >= required_start_pos)
                        && (current_ref_pos <= required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                        let (p1, p2) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, required_start_pos, THREE_BASE_CONTEXT_READ_LENGTH);
                        temp_read_vec = [temp_read_vec, readunwrapped.seq().as_bytes()[p1..p2].to_vec()].concat();
                    }
                    current_ref_pos += temp_int;
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'H' => {
                    temp_character_vec = vec![];
                },
                'S' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'I' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'N' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    //println!("N RUNNIN at pos {} for {}", current_ref_pos, temp_int); 
                    if (current_ref_pos + temp_int >= required_start_pos)
                        && (current_ref_pos <= required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                        let (p1, p2) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, required_start_pos, THREE_BASE_CONTEXT_READ_LENGTH);
                        for _ in p1..p2 {
                            temp_read_vec.push('X' as u8);
                        }
                    }
                    current_ref_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'D' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    //println!("D RUNNIN at pos {} for {}", current_ref_pos, temp_int); 
                    if (current_ref_pos + temp_int >= required_start_pos)
                        && (current_ref_pos <= required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                        let (p1, p2) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, required_start_pos, THREE_BASE_CONTEXT_READ_LENGTH);
                        for _ in p1..p2 {
                            temp_read_vec.push('X' as u8);
                        }
                    }
                    current_ref_pos += temp_int;
                    temp_character_vec = vec![];
                },
                _ => {
                    temp_character_vec.push(*character as char);
                },
            }
        }
        // take care of cases when string is not long enough to cover the whole ref required
        read_vec.push(String::from_utf8(temp_read_vec).expect("Found invalid UTF-8"));
    }
    read_vec
}

fn read_fai_file (start_pos: usize, chromosone: &String, reader: &mut faidx::Reader) -> String {
    //println!("The ref sequenece");
    //println!("{}", fasta_reader.fetch_seq_string(chromosone, start_pos, start_pos + THREE_BASE_CONTEXT_READ_LENGTH).unwrap());
    reader.fetch_seq_string(chromosone, start_pos, start_pos + THREE_BASE_CONTEXT_READ_LENGTH).unwrap()
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