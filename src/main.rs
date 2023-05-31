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
use rust_htslib::bam::{Read as BamRead, IndexedReader};
use rust_htslib::bcf::{Reader, Read as BcfRead};
use std::convert::TryFrom;

const SEED: u64 = 2;
const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = 0;
const MATCH: i32 = 2;
const MISMATCH: i32 = -2;
const RANDOM_SEQUENCE_LENGTH: usize = 1000;
const NUMBER_OF_RANDOM_SEQUENCES: usize = 5;
const THREE_BASE_CONTEXT_READ_LENGTH: usize = 1000;

fn main() {
    //homopolymer_test();
    //let seqvec = get_random_sequences_from_generator(RANDOM_SEQUENCE_LENGTH, NUMBER_OF_RANDOM_SEQUENCES, SEED);
    //pipeline_quality_score();
    //pipeline_3base_context();
    read_vcf_file();
}

fn pipeline_quality_score_error_graph () {
    // get the errors from himut vcf (no somatic mutations in this file)

    // get the quality scores in the ccs
    
    // get the actual counts of the errors

    // draw graph

}

fn read_vcf_file() {
    let path = &"data/output.vcf.gz";
    let mut bcf = Reader::from_path(path).expect("Error opening file.");
    // iterate through each row of the vcf body.
    for (_, record_result) in bcf.records().enumerate() {
        let record = record_result.expect("Fail to read record");
        let mut s = String::new();
        println!("{:?}", record.desc());
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
        println!("sample count: {}", n_ref.len());
        println!("num of referce allele: {}\nnum of alt allele: {}\nnum of missing allele: {}\n\n", n_ref[0], n_alt[0], n_missing[0]);
    }
}

fn pipeline_quality_score () {
    // find the error position
    let error_pos = 1000000;
    let error_chr = String::from("chr1");
    // get the reads corrosponding to the position
    let reads = get_read_and_readnames_from_bam(&error_pos, &error_chr);
    // get the sub reads of the reads // do poa with the reads
    for read in &reads {
        get_subreads_from_readname(&read.1, &error_pos, &error_chr);

    }
    // check the quality scores at that location or if fixed
}

fn get_read_and_readnames_from_bam (error_pos: &usize, error_chr: &String) -> Vec<(String, String, usize)> {
    let mut reads_names_and_errorpos: Vec<(String, String, usize)> = vec![];
    let path = &"data/merged.bam";
    let mut bam_reader = IndexedReader::from_path(path).unwrap();
    bam_reader.fetch((error_chr, *error_pos as i64, *error_pos as i64 + 1)).unwrap();
    for read in bam_reader.records() {
        let readunwrapped = read.unwrap();
        let read_index = error_pos - readunwrapped.pos() as usize;
        let read_name = String::from_utf8(readunwrapped.qname().to_vec()).expect("");
        let read_string = String::from_utf8(readunwrapped.seq().as_bytes().to_vec()).expect("");
        reads_names_and_errorpos.push((read_string.clone(),read_name.clone(), read_index));
    }
    reads_names_and_errorpos
}

fn get_subreads_from_readname (read_name: &String, error_pos: &usize, error_chr: &String) -> Vec<String>  {
    let mut subreads = vec![];
    // open the bam file corrosponding to the read_name
    println!("{}", read_name);
    let mut split_text_iter = (read_name.split("/")).into_iter();
    let chip_name = split_text_iter.next().unwrap();
    let consensus_name = split_text_iter.next().unwrap();
    // get the subreads corrosponding to read_name and position
    let bam_file_path = format!("/Users/wmw0016/Documents/mount5/data1/hifi_consensus/try2/{}.subreads.mapped.bam", chip_name);
    let mut bam_reader = IndexedReader::from_path(bam_file_path).unwrap();
    bam_reader.fetch((error_chr, *error_pos as i64, *error_pos as i64 + 1)).unwrap();
    for read in bam_reader.records() {
        let readunwrapped = read.unwrap();
        let read_name = String::from_utf8(readunwrapped.qname().to_vec()).expect("");
        let mut read_name_iter = read_name.split("/");
        let subread_chip = read_name_iter.next().unwrap();
        let subread_name = read_name_iter.next().unwrap();
        if subread_name == consensus_name {
            subreads.push(String::from_utf8(readunwrapped.seq().as_bytes().to_vec()).expect(""));
            println!("{} {}", subread_chip, subread_name);
        }
    }
    subreads
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
            let mut bam_reader = IndexedReader::from_path(path).unwrap();
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



fn read_bam_file (required_start_pos: usize, chromosone: &String, reader: &mut IndexedReader) -> Vec<String> {
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
                        let p1;
                        let p2;
                        // case when both end partial match is required
                        if (current_ref_pos < required_start_pos) && (current_ref_pos + temp_int > required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                            p1 = required_start_pos - current_ref_pos + current_read_pos;
                            p2 = current_read_pos + required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH - current_ref_pos;
                        }
                        // case when start partial matches are required
                        else if (current_ref_pos < required_start_pos) && (current_ref_pos + temp_int >= required_start_pos) {
                            p1 = required_start_pos - current_ref_pos + current_read_pos;
                            p2 = current_read_pos + temp_int;
                        }
                        // case when end partial matches are required
                        else if (current_ref_pos <= required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) && (current_ref_pos + temp_int > required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                            p1 = current_read_pos;
                            p2 = current_read_pos + required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH - current_ref_pos;
                        }
                        // normal case
                        else {
                            p1 = current_read_pos;
                            p2 = current_read_pos + temp_int;
                        }
                        temp_read_vec = [temp_read_vec, readunwrapped.seq().as_bytes()[p1..p2].to_vec()].concat();
                    }
                    current_ref_pos += temp_int;
                    current_read_pos += temp_int;
                    temp_character_vec = vec![];
                },
                'H' => {
                    //let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    //let temp_int = temp_string.parse::<usize>().unwrap();
                    //current_read_pos += temp_int;
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
                        let p1;
                        let p2;
                        // case when both end partial match is required
                        if (current_ref_pos < required_start_pos) && (current_ref_pos + temp_int > required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                            p1 = required_start_pos - current_ref_pos + current_read_pos;
                            p2 = current_read_pos + required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH - current_ref_pos;
                        }
                        // case when start partial matches are required
                        else if (current_ref_pos < required_start_pos) && (current_ref_pos + temp_int >= required_start_pos) {
                            p1 = required_start_pos - current_ref_pos + current_read_pos;
                            p2 = current_read_pos + temp_int;
                        }
                        // case when end partial matches are required
                        else if (current_ref_pos <= required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) && (current_ref_pos + temp_int > required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                            p1 = current_read_pos;
                            p2 = current_read_pos + required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH - current_ref_pos;
                        }
                        // normal case
                        else {
                            p1 = current_read_pos;
                            p2 = current_read_pos + temp_int;
                        }
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
                        let p1;
                        let p2;
                        // case when both end partial match is required
                        if (current_ref_pos < required_start_pos) && (current_ref_pos + temp_int > required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                            p1 = required_start_pos - current_ref_pos + current_read_pos;
                            p2 = current_read_pos + required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH - current_ref_pos;
                        }
                        // case when start partial matches are required
                        else if (current_ref_pos < required_start_pos) && (current_ref_pos + temp_int >= required_start_pos) {
                            p1 = required_start_pos - current_ref_pos + current_read_pos;
                            p2 = current_read_pos + temp_int;
                        }
                        // case when end partial matches are required
                        else if (current_ref_pos <= required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) && (current_ref_pos + temp_int > required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH) {
                            p1 = current_read_pos;
                            p2 = current_read_pos + required_start_pos + THREE_BASE_CONTEXT_READ_LENGTH - current_ref_pos;
                        }
                        // normal case
                        else {
                            p1 = current_read_pos;
                            p2 = current_read_pos + temp_int;
                        }
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
    for _read in &read_vec {
        //println!("{}",read);
        //println!("read ln {}", read.len());
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