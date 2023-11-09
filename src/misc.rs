extern crate bio;
use bio::alignment::pairwise::banded::Aligner as BandedDP;
use libm::sqrt;

use crate::alignment::poabandedsmarter::Aligner;
use crate::quality::topology_cut::get_consensus_parallel_bases;
use petgraph::{Graph, Directed, graph::NodeIndex};
use petgraph::dot::Dot;
use petgraph::visit::Topo;
use petgraph::Direction::Outgoing;
use rust_htslib::bam::{Read as BamRead, IndexedReader as BamIndexedReader, record::Aux};
use rust_htslib::bcf::{Reader, Read as BcfRead};
use rust_htslib::faidx;
use std::{fs::OpenOptions, io::{prelude::*}};
use std::io::BufReader;
use std::fs::File;
use std::fs::remove_file;
use std::io::SeekFrom;
use std::fs::read_dir;
use std::fs::create_dir_all;
use std::{io};
use std::fs::read_to_string;

const SEED: u64 = 3;
const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = 0;
const MATCH: i32 = 2;
const MISMATCH: i32 = -2;
const RANDOM_SEQUENCE_LENGTH: usize = 6;
const NUMBER_OF_RANDOM_SEQUENCES: usize = 20;
const THREE_BASE_CONTEXT_READ_LENGTH: usize = 2;
const NUM_OF_ITER_FOR_ZOOMED_GRAPHS: usize = 4;
const DATA_PATH: &str = "/data1/hifi_consensus/try2/";
const READ_BAM_PATH: &str = "/data1/hifi_consensus/try2/merged.bam";
const INTERMEDIATE_PATH: &str = "/data1/hifi_consensus/intermediate";
const CONFIDENT_PATH: &str = "/data1/GiaB_benchmark/HG001_GRCh38_1_22_v4.2.1_benchmark.bed";
const REF_GENOME_PATH: &str = "/data1/GiaB_benchmark/GRCh38.fa";
const RESULT_WRITE_PATH: &str = "/data1/hifi_consensus/all_data/chr2_7base_data";
const DEEPVARIANT_PATH: &str = "/data1/hifi_consensus/try3/hg38.PD47269d.minimap2_ccs.deepvariant_1.1.0.vcf";
const WRONG_ERROR_FILE_PATH: &str = "/data1/hifi_consensus/all_data/chr2_errors.txt";
const HIMUT_PATH: &str = "/data1/hifi_consensus/try3/test.vcf";
const BAND_SIZE: i32 = 100;
const MAX_NODES_IN_POA: usize = 75_000;
const SKIP_SCORE: i32 = 6_000;

pub fn pipeline_load_graph_get_topological_parallel_bases (chromosone: &str, start: usize, end: usize, thread_id: usize) {
    let mut index_thread = 0;
    let mut skip_thousand = false;
    let mut skip_index = 0;
    for process_location in start..end {
        // skip thousand when same found
        if skip_thousand {
            skip_index += 1;
            if skip_index > 5000 {
                skip_thousand = false;
                skip_index = 0;
            }
            else {
                continue;
            }
        }
        println!("Thread {}: Chr {} Loc {}, tasks_done {} NEW LOCATION", thread_id, chromosone, process_location, index_thread);
        // get the string and the name
        let seq_name_qual_and_errorpos_vec = get_corrosponding_seq_name_location_quality_from_bam(process_location, &chromosone.to_string(), &'X');
        let mut all_skipped = true;
        for seq_name_qual_and_errorpos in &seq_name_qual_and_errorpos_vec {
            // check if the css file is already available
            let check_file = format!("{}_parallel.txt", &seq_name_qual_and_errorpos.1);
            if check_file_availability(&check_file, INTERMEDIATE_PATH) {
                //println!("Thread {}: Required CSS File Available, skipping..", thread_id);
                continue;
            }
            // check if graph is available, if available load all the data
            let check_file = format!("{}_graph.txt", &seq_name_qual_and_errorpos.1);
            if check_file_availability(&check_file, INTERMEDIATE_PATH) {
                //println!("Thread {}: Required File not Available, Graph Available, processing..", thread_id);
            }
            // if both not available
            else {
                //println!("Thread {}: Nothing is available, continuing..", thread_id);
                continue;
            }
            all_skipped = false;
            // find the subreads of that ccs
            let sub_reads = get_the_subreads_by_name_sam(&seq_name_qual_and_errorpos.1);
            // skip if no subreads, errors and stuff
            if sub_reads.len() == 0 {
                continue;
            }
            let calculated_graph = load_the_graph(check_file);
            let (calculated_consensus, calculated_topology) = get_consensus_from_graph(&calculated_graph); //just poa
            let parallel_bases_vec = get_consensus_parallel_bases(sub_reads.len(), &calculated_consensus, &calculated_topology, &calculated_graph, thread_id);
            // match the calculated consensus to the original consensus and get the required indices
            let calc_cons_id = get_redone_consensus_matched_positions(&seq_name_qual_and_errorpos.0, &calculated_consensus);
            for (index, pacbio_base) in seq_name_qual_and_errorpos.0.as_bytes().to_vec().iter().enumerate() {
                let mut pacbio_str = format!("OK({})", *pacbio_base as char);
                let parallel_bases;
                if calc_cons_id[index].1 == 0 {
                    parallel_bases = vec![1, 1, 1, 1]; //deletion
                    pacbio_str = "DEL()".to_string();
                }
                else if calc_cons_id[index].1 == 1 {
                    parallel_bases = parallel_bases_vec[calc_cons_id[index].0].clone(); //normal
                }
                else {
                    parallel_bases = parallel_bases_vec[calc_cons_id[index].0].clone(); //subsitution the value corrospond to the sub
                    pacbio_str = format!{"SB({})", calc_cons_id[index].1 as char};
                }
                let write_string = format!("{} {:?}\n", pacbio_str, parallel_bases);
                let write_file = format!("{}/{}_parallel.txt", INTERMEDIATE_PATH, &seq_name_qual_and_errorpos.1);
                write_string_to_file(&write_file, &write_string);
            } 
            index_thread += 1;
            println!("Thread {}: Chr {} Loc {}, tasks_done {}", thread_id, chromosone, process_location, index_thread);
        }
        if all_skipped {
            skip_thousand = true;
        }
    }
}

fn get_redone_consensus_matched_positions (pacbio_consensus: &String, calculated_consensus: &Vec<u8>) -> Vec<(usize, u8)> {
    let mut consensus_matched_indices: Vec<(usize, u8)> = vec![];
    let pacbio_consensus_vec: Vec<u8> = pacbio_consensus.bytes().collect();
    let score_func = |a: u8, b: u8| if a == b { 2i32 } else { -2i32 };
    let k = 8; // kmer match length
    let w = 20; // Window size for creating the band
    let mut aligner = BandedDP::new(-2, -2, score_func, k, w);
    let alignment = aligner.global(&calculated_consensus, &pacbio_consensus_vec);
    let temp_score = alignment.score;
    //let (alignment, temp_score) = pairwise(&calculated_consensus, &pacbio_consensus_vec, MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND, 0);
    let mut calc_index = 0;
    for op in alignment.operations {
        match op {
            bio::alignment::AlignmentOperation::Match => {
                consensus_matched_indices.push((calc_index, 1));
                calc_index += 1;
            },
            bio::alignment::AlignmentOperation::Subst => {
                let base = calculated_consensus[calc_index];
                consensus_matched_indices.push((calc_index, base));
                calc_index += 1;
            },
            bio::alignment::AlignmentOperation::Del => {
                consensus_matched_indices.push((calc_index, 0));
            },
            bio::alignment::AlignmentOperation::Ins => {
                calc_index += 1;
            },
            _ => {},
        }
    }
    println!("score {}", temp_score);
    consensus_matched_indices
}

pub fn get_all_data_for_ml (chromosone: &str, start: usize, end: usize, thread_id: usize) {
    let mut position_base = start;
    'bigloop: loop {
        if position_base % 1000 == 0 {
            println!("Thread ID: {} Position {}", thread_id, position_base);
        }
        let seq_name_qual_and_errorpos_vec = get_corrosponding_seq_name_location_quality_from_bam(position_base, &chromosone.to_string(), &'X');
        // get the three base context
        let mut ref_sevenbase_context = "".to_string();
        if seq_name_qual_and_errorpos_vec.len() > 0 {
            let mut fai_reader = faidx::Reader::from_path(REF_GENOME_PATH).unwrap();
            ref_sevenbase_context = read_fai_get_ref_context(position_base - 3,7,  &chromosone.to_string(), &mut fai_reader);
        }
        for seq_name_qual_and_errorpos in &seq_name_qual_and_errorpos_vec {
            let char_sequence: Vec<char> = seq_name_qual_and_errorpos.0.chars().collect::<Vec<_>>();
            let mut char_7base_context: Vec<char> = vec![];
            // when 7 base context above 0 this is very dumb re write
            if seq_name_qual_and_errorpos.3 >= 3 {
                char_7base_context.push(char_sequence[seq_name_qual_and_errorpos.3 - 3]);
            }
            else {
                char_7base_context.push('X');
            }
            if seq_name_qual_and_errorpos.3 >= 2 {
                char_7base_context.push(char_sequence[seq_name_qual_and_errorpos.3 - 2]);
            }
            else{
                char_7base_context.push('X');
            }
            if seq_name_qual_and_errorpos.3 >= 1 {
                char_7base_context.push(char_sequence[seq_name_qual_and_errorpos.3 - 1]);
            }
            else {
                char_7base_context.push('X');
            }
            char_7base_context.push(char_sequence[seq_name_qual_and_errorpos.3]);
            // when len is greater than 7 base context
            if seq_name_qual_and_errorpos.0.len() > (1 + seq_name_qual_and_errorpos.3) {
                char_7base_context.push(char_sequence[seq_name_qual_and_errorpos.3 + 1]);
            }
            // when len is less than 7 base context
            else {
                char_7base_context.push('X');
            }
            if seq_name_qual_and_errorpos.0.len() > (2 + seq_name_qual_and_errorpos.3) {
                char_7base_context.push(char_sequence[seq_name_qual_and_errorpos.3 + 2]);
            }
            else{
                char_7base_context.push('X');
            }
            if seq_name_qual_and_errorpos.0.len() > (3 + seq_name_qual_and_errorpos.3) {
                char_7base_context.push(char_sequence[seq_name_qual_and_errorpos.3 + 3]);
            }
            else{
                char_7base_context.push('X');
            }
            let read_sevenbase_context = char_7base_context.iter().collect::<String>();
            let quality = seq_name_qual_and_errorpos.2;
            // error is here
            let parallel_stuff;
            // check if the file is already available
            let file_name = format!("{}{}", seq_name_qual_and_errorpos.1, "_parallel.txt");
            if check_file_availability(&file_name, INTERMEDIATE_PATH) {
                let available_file_path = format!("{}/{}", INTERMEDIATE_PATH, file_name);
                parallel_stuff = get_parallel_bases_from_file(&available_file_path, seq_name_qual_and_errorpos.3);
            }
            else {
                continue;
            }
            let read_position = seq_name_qual_and_errorpos.3;
            let read_len =  seq_name_qual_and_errorpos.0.len();
            // write data
            let write_string = format!("{} {} {} : {} {} {} {} {:?}", position_base, ref_sevenbase_context, quality, read_position, read_len, read_sevenbase_context, parallel_stuff, seq_name_qual_and_errorpos.4);
            println!("{}", write_string);
            let write_file = format!("{}/{}_mldata.txt", RESULT_WRITE_PATH, thread_id);
            //write_string_to_file(&write_file, &write_string);
        }
        position_base += 1;
        if position_base > end {
            break 'bigloop;
        }
    }
}

pub fn create_himut_list () {
    // get the error locations
    let error_locations = get_himut_info_from_himut_vcf (); //chromosone, location, ref allele, alt allele
    for error_location in error_locations {
        let error_string = format!("{} {} {} -> {}\n", error_location.0, error_location.1, error_location.2, error_location.3);
        let write_file = format!("/data1/hifi_consensus/all_data/filters/himut_data.txt");
        write_string_to_file(&write_file, &error_string);
    }
}

pub fn get_himut_info_from_himut_vcf () -> Vec<(String, usize, String, String)> {
    let mut error_locus_vec: Vec<(String, usize, String, String)> = vec![]; //chromosone, position, ref, alt
    let mut bcf = Reader::from_path(HIMUT_PATH).expect("Error opening file.");
    // iterate through each row of the vcf body.
    for (_, record_result) in bcf.records().enumerate() {
        let record = record_result.expect("Fail to read record");
        let mut allele_vec: Vec<String> = vec![];
        // only pass filters are accepted
        //if record.has_filter("PASS".as_bytes()) == true {
            for allele in record.alleles() {
                allele_vec.push(std::str::from_utf8(allele).unwrap().to_owned());
            }
            let temp_str = record.desc();
            let mut split_text_iter = (temp_str.split(":")).into_iter();
            let chromosone = format!("{}", split_text_iter.next().unwrap());
            error_locus_vec.push((chromosone.to_string(), record.pos() as usize, allele_vec[0].clone(), allele_vec[1].clone()));
            println!("{} {} {} -> {}", chromosone.to_string(), record.pos() as usize, allele_vec[0], allele_vec[1]);
        //}
    }
    println!("number of errors = {}", error_locus_vec.len());
    error_locus_vec
}

pub fn create_depth_indel_list (chromosone: &str, start: usize, end: usize, thread_id: usize) {
    // go though the locations
    for position_base in start..end {
        if position_base % 1000 == 0 {
            println!("Position {}", position_base);
        }
        let cause_value = get_info_from_bam(position_base, &chromosone.to_string());
        if cause_value.0 == 1 {
            let write_string = format!("{} {} {}\n", chromosone, position_base, cause_value.1);
            let write_file = format!("{}/{}_depth.txt", RESULT_WRITE_PATH, thread_id);
            write_string_to_file(&write_file, &write_string);
        }
        else if cause_value.0 == 2 {
            let write_string = format!("{} {} {}/{}\n", chromosone, position_base, cause_value.2, cause_value.1);
            let write_file = format!("{}/{}_indel.txt", RESULT_WRITE_PATH, thread_id);
            write_string_to_file(&write_file, &write_string);
        }
    }
}

fn get_info_from_bam (error_pos: usize, error_chr: &String) -> (usize, usize, usize) {
    let path = &READ_BAM_PATH;
    let mut bam_reader = BamIndexedReader::from_path(path).unwrap();
    match bam_reader.fetch((error_chr, error_pos as i64, error_pos as i64 + 1)) {
        Ok(_x) => {},
        _ => {return (0, 0, 0);},
    }
    let mut overlap_indel_count: usize = 0;
    let mut depth_count: usize = 0;
    'read_loop: for read in bam_reader.records() {
        depth_count += 1;
        let readunwrapped = read.unwrap();
        // decode the cigar string
        // get the location from the cigar processing
        let mut temp_character_vec: Vec<char> = vec![];
        // get the read start position
        let read_start_pos = readunwrapped.pos() as usize;
        let mut current_ref_pos = read_start_pos;
        let mut _current_read_pos = 0;
        let mut last_one_del = false;
        for character in readunwrapped.cigar().to_string().as_bytes() {
            match *character as char {
                'M' => {         
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    if current_ref_pos + 1 == error_pos {
                        if last_one_del {
                            overlap_indel_count += 1;
                        }
                    }
                    if (current_ref_pos + temp_int >= error_pos)
                        && (current_ref_pos <= error_pos + 1) {
                        continue 'read_loop;
                    }
                    current_ref_pos += temp_int;
                    _current_read_pos += temp_int;
                    temp_character_vec = vec![];
                    last_one_del = false;
                },
                'H' => {
                    temp_character_vec = vec![];
                    last_one_del = false;
                },
                'S' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    _current_read_pos += temp_int;
                    temp_character_vec = vec![];
                    last_one_del = false;
                },
                'I' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    _current_read_pos += temp_int;
                    temp_character_vec = vec![];
                    last_one_del = true;
                },
                'N' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    if (current_ref_pos + temp_int >= error_pos)
                        && (current_ref_pos <= error_pos + 1) {
                        continue 'read_loop;
                    }
                    current_ref_pos += temp_int;
                    temp_character_vec = vec![];
                    last_one_del = false;
                },
                'D' => {
                    let temp_string: String = temp_character_vec.clone().into_iter().collect();
                    let temp_int = temp_string.parse::<usize>().unwrap();
                    if (current_ref_pos + temp_int >= error_pos)
                        && (current_ref_pos <= error_pos + 1) {
                        continue 'read_loop;
                    }
                    current_ref_pos += temp_int;
                    temp_character_vec = vec![];
                    last_one_del = false;
                },
                _ => {
                    temp_character_vec.push(*character as char);
                },
            }
        }
    }
    drop(bam_reader);
    let mut cause: usize = 0;
    if depth_count as f64 > (30.0 + (4.0 * sqrt(30.0))) {
        cause = 1;
    }
    else if (overlap_indel_count < depth_count) && (overlap_indel_count > 0) {
        cause = 2;
    }
    (cause, depth_count, overlap_indel_count)
}

pub fn create_confidence_list () {
    // get the confident locations
    let confident_locations = get_confident_locations_from_file (); //chromosone, location, ref allele, alt allele
    for confident_location in confident_locations {
        let error_string = format!("{} {} ~ {}\n", confident_location.0, confident_location.1, confident_location.2);
        let write_file = format!("/data1/hifi_consensus/all_data/filters/confident_data.txt");
        write_string_to_file(&write_file, &error_string);
    }
}

fn get_confident_locations_from_file () -> Vec<(String, usize, usize)> {
    let mut location_vec: Vec<(String, usize, usize)> = vec![];
    let file_path = CONFIDENT_PATH;
    let f = File::open(&file_path).unwrap();
    let mut reader = BufReader::new(f);
    let mut buffer = String::new();
    loop {
        buffer.clear();
        match reader.read_line(&mut buffer) {
            Ok(_) => {
                let mut split_text_iter = (buffer.split("\t")).into_iter();
                let chromosone ; 
                match split_text_iter.next() {
                    Some(x) => {chromosone = x.to_string();},
                    None => {break;},
                }
                let start_string;
                match split_text_iter.next() {
                    Some(x) => {start_string = x.to_string();},
                    None => {break;},
                }
                let mut end_string;
                match split_text_iter.next() {
                    Some(x) => {end_string = x.to_string();},
                    None => {break;},
                };
                end_string.pop();
                let end = end_string.parse::<usize>().unwrap();
                let start = start_string.parse::<usize>().unwrap();
                //println!("{} {} {}",chromosone, start, end);
                location_vec.push((chromosone, start, end));
            },
            Err(_) => {break;},
        };
    }
    location_vec
}

pub fn concancate_files () {
    let mut output = File::create("/data1/hifi_consensus/all_data/chr2.txt").unwrap();
    //let paths = read_dir("data/chr21/").unwrap();
    let parent_path = "/data1/hifi_consensus/all_data/chr2_data/";
    let mut path_array = vec![];
    for i in 0..70 {
        path_array.push(format!("{}{}_mldata.txt",parent_path, i))
    }
    // convert path to string and save in array
    for path in path_array {
        let input = File::open(path.clone());
        match input {
            Ok(mut x) => {
                io::copy(&mut x, &mut output).unwrap();
                println!("{}", path);
            },
            _ => {},
        }
    }
    println!("done");
}

pub fn create_germline_list () {
    // get the error locations
    let error_locations = get_germline_info_from_deepvariant_vcf (); //chromosone, location, ref allele, alt allele
    for error_location in error_locations {
        let error_string = format!("{} {} {} -> {}\n", error_location.0, error_location.1, error_location.2, error_location.3);
        let write_file = format!("/data1/hifi_consensus/all_data/filters/germline_data.txt");
        write_string_to_file(&write_file, &error_string);
    }
}

pub fn get_germline_info_from_deepvariant_vcf () -> Vec<(String, usize, String, String)> {
    let mut error_locus_vec: Vec<(String, usize, String, String)> = vec![]; //chromosone, position, ref, alt
    let mut bcf = Reader::from_path(DEEPVARIANT_PATH).expect("Error opening file.");
    // iterate through each row of the vcf body.
    for (_, record_result) in bcf.records().enumerate() {
        let record = record_result.expect("Fail to read record");
        let mut allele_vec: Vec<String> = vec![];
        // only pass filters are accepted
        if record.has_filter("PASS".as_bytes()) == true {
            for allele in record.alleles() {
                allele_vec.push(std::str::from_utf8(allele).unwrap().to_owned());
            }
            let temp_str = record.desc();
            let mut split_text_iter = (temp_str.split(":")).into_iter();
            let chromosone = format!("{}", split_text_iter.next().unwrap());
            error_locus_vec.push((chromosone.to_string(), record.pos() as usize, allele_vec[0].clone(), allele_vec[1].clone()));
            println!("{} {} {} -> {}", chromosone.to_string(), record.pos() as usize, allele_vec[0], allele_vec[1]);
        }
    }
    println!("number of errors = {}", error_locus_vec.len());
    error_locus_vec
}

fn get_redone_consensus_error_position (pacbio_consensus: &String, calculated_consensus: &Vec<u8>, pacbio_error_position: usize) -> usize {
    let mut consensus_match_invalid_indices: Vec<usize> = vec![];
    let mut aligned_pacbio_scores_vec: Vec<usize> = vec![];
    let mut aligned_pacbio_bases:Vec<u8> = vec![];
    let pacbio_consensus_vec: Vec<u8> = pacbio_consensus.bytes().collect();
    let score_func = |a: u8, b: u8| if a == b { 2i32 } else { -2i32 };
    let k = 8; // kmer match length
    let w = 20; // Window size for creating the band
    let mut aligner = BandedDP::new(-2, -2, score_func, k, w);
    let alignment = aligner.global(&calculated_consensus, &pacbio_consensus_vec);
    let mut pacbio_index = 0;
    let mut calc_index = 0;
    let mut calc_error_position: usize = 0;
    for op in alignment.operations {
        match op {
            bio::alignment::AlignmentOperation::Match => {
                pacbio_index += 1;
                calc_index += 1;
            },
            bio::alignment::AlignmentOperation::Subst => {
                aligned_pacbio_bases.push(pacbio_consensus_vec[pacbio_index]);
                consensus_match_invalid_indices.push(calc_index);
                pacbio_index += 1;
                calc_index += 1;
            },
            bio::alignment::AlignmentOperation::Del => {
                pacbio_index += 1;
            },
            bio::alignment::AlignmentOperation::Ins => {
                consensus_match_invalid_indices.push(calc_index);
                aligned_pacbio_bases.push(126);
                aligned_pacbio_scores_vec.push(33);
                calc_index += 1;
            },
            _ => {},
        }
        if pacbio_error_position == pacbio_index {
            calc_error_position = calc_index
        }
    }
    println!("Pacbio error position {} corrosponds to calculated error position {}", pacbio_error_position, calc_error_position);
    calc_error_position
}

pub fn pipeline_save_the_graphs (chromosone: &str, start: usize, end: usize, thread_id: usize) {
    let mut big_file_skip_count = 0;
    let mut index_thread = 0;
    let mut skip_thousand = false;
    let mut skip_index = 0;
    'bigloop: for process_location in start..end {
        // skip thousand when same found
        if skip_thousand {
            skip_index += 1;
            if skip_index > 5000 {
                skip_thousand = false;
                skip_index = 0;
            }
            else {
                continue;
            }
        }
        println!("NEW LOCATION, Thread {}: Chr {} Loc {}, tasks_done {} skipped {}", thread_id, chromosone, process_location, index_thread, big_file_skip_count);
        // get the string and the name
        let seq_name_qual_and_errorpos_vec = get_corrosponding_seq_name_location_quality_from_bam(process_location, &chromosone.to_string(), &'X');
        let mut all_skipped = true;
        for seq_name_qual_and_errorpos in &seq_name_qual_and_errorpos_vec {
            println!("Thread {}: Chr {} Loc {} Processing ccs file: {}", thread_id, chromosone, process_location, seq_name_qual_and_errorpos.1);
            // check if the graph is already available
            let check_file = format!("{}_graph.txt", &seq_name_qual_and_errorpos.1);
            if check_file_availability(&check_file, INTERMEDIATE_PATH) {
                println!("Thread {}: File Available, skipping", thread_id);
                continue;
            }
            // if not available do poa and make a file
            // find the subreads of that ccs
            let mut sub_reads = get_the_subreads_by_name_sam(&seq_name_qual_and_errorpos.1);
            // skip if no subreads, errors and stuff
            if sub_reads.len() == 0 {
                continue;
            }
            all_skipped = false;
            // filter out the long reads and rearrange the reads
            sub_reads = reverse_complement_filter_and_rearrange_subreads(&sub_reads);
            // reverse if score is too low
            sub_reads = check_the_scores_and_change_alignment(sub_reads, &seq_name_qual_and_errorpos.0);
            if sub_reads.len() == 0 {
                skip_thousand = true;
                continue 'bigloop;
            }

            sub_reads.insert(0, seq_name_qual_and_errorpos.0.clone());
            // do poa with the read and subreads, get the poa and consensus
            let mut sequence_number: usize = 0;
            let mut aligner = Aligner::new(MATCH, MISMATCH, GAP_OPEN, &sub_reads[0].as_bytes().to_vec(), BAND_SIZE);
            
            for sub_read in &sub_reads {
                if sequence_number != 0 {
                    aligner.global(&sub_read.as_bytes().to_vec()).add_to_graph();
                }
                let node_num = aligner.graph().node_count();
                if node_num > MAX_NODES_IN_POA {
                    println!("NUM OF NODES {} TOO BIG, SKIPPING TOTAL SKIPPED: {} ", node_num, big_file_skip_count + 1);
                    big_file_skip_count += 1;
                    skip_thousand = true;
                    continue 'bigloop;
                }
                sequence_number += 1;
                println!("Thread {}: Sequence {} processed", thread_id, sequence_number);
            }
            let calculated_graph: &Graph<u8, i32, Directed, usize> = aligner.graph();
            save_the_graph(calculated_graph, &seq_name_qual_and_errorpos.1);
            index_thread += 1;
        }
        if all_skipped {
            skip_thousand = true;
        }
    }
}

fn get_consensus_from_graph(graph: &Graph<u8, i32, Directed, usize>) -> (Vec<u8>, Vec<usize>) {
    let mut output: Vec<u8> = vec![];
    let mut topopos: Vec<usize> = vec![];
    let mut topo = Topo::new(graph);
    let mut topo_indices = Vec::new();
    let mut max_index = 0;
    let mut max_score = 0.0;

    while let Some(node) = topo.next(graph) {
        topo_indices.push(node);
        if max_index < node.index(){
            max_index = node.index();
        }
    }
    topo_indices.reverse();
    //define score and nextinpath vectors with capacity of num nodes.
    let mut weight_scores: Vec<i32> = vec![0; max_index + 1];
    let mut scores: Vec<f64> = vec![0.0; max_index + 1];
    let mut next_in_path: Vec<usize> = vec![0; max_index + 1];
    //iterate thorugh the nodes in reverse
    for node in topo_indices{
        let mut best_weight_score_edge: (i32, f64, usize) = (-1 , -1.0, 123456789);
        let mut neighbour_nodes = graph.neighbors_directed(node, Outgoing);
        while let Some(neighbour_node) = neighbour_nodes.next() {
            let mut edges = graph.edges_connecting(node, neighbour_node);
            let mut weight: i32 = 0;
            while let Some(edge) = edges.next() {
                weight += edge.weight().clone();
            }
            let weight_score_edge = (weight, scores[neighbour_node.index()], neighbour_node.index());
            if weight_score_edge > best_weight_score_edge{
                best_weight_score_edge = weight_score_edge;
            }
        }
        //save score and traceback
        if best_weight_score_edge.0 as f64 + best_weight_score_edge.1 > max_score{
            max_score = best_weight_score_edge.0 as f64 + best_weight_score_edge.1;
        }
        scores[node.index()] = best_weight_score_edge.0 as f64 + best_weight_score_edge.1;
        next_in_path[node.index()] = best_weight_score_edge.2;
        weight_scores[node.index()] = best_weight_score_edge.0;
    }
    let mut pos = scores.iter().position(|&r| r == max_score).unwrap();
    //calculate the start weight score
    let mut consensus_started: bool = false;
    let weight_average = scores[pos] / scores.len() as f64;
    let weight_threshold = weight_average as i32 / 2; 
    while pos != 123456789 {
        //continue if starting weight score is too low
        if consensus_started == false && weight_scores[pos] < weight_threshold {
            pos = next_in_path[pos];
            continue;
        }
        //println!("current {} {}", pos, next_in_path[pos]);
        consensus_started = true;
        topopos.push(pos as usize);
        output.push(graph.raw_nodes()[pos].weight);
        pos = next_in_path[pos];
    }
    (output, topopos)
}

pub fn write_string_to_newfile (file_name: &String, input_string: &String) {
    let path = std::path::Path::new(&file_name);
    let prefix = path.parent().unwrap();
    match remove_file(path) {
        Ok(_) => {},
        Err(_) => {}
    };
    create_dir_all(prefix).unwrap();
    let mut file = OpenOptions::new().create(true).write(true).open(file_name).unwrap();
    writeln!(file, "{}", input_string).expect("result file cannot be written");
}

fn save_the_graph (graph: &Graph<u8, i32, Directed, usize>, file_name: &String) {
    let mut write_string = "".to_string();
    let write_path = format!("{}/{}_graph.txt", INTERMEDIATE_PATH, file_name);
    let mut node_iterator = graph.node_indices();
    while let Some(node) = node_iterator.next() {
        let node_index = node.index();
        let base = graph.raw_nodes()[node_index].weight;
        let mut neighbours:Vec<(usize, i32)> = vec![];
        let mut neighbour_nodes = graph.neighbors_directed(node, Outgoing);
        let mut neighbour_string = "".to_string();
        while let Some(neighbour_node) = neighbour_nodes.next() {
            let mut edges = graph.edges_connecting(node, neighbour_node);
            let mut weight: i32 = 0;
            while let Some(edge) = edges.next() {
                weight += edge.weight().clone();
            }
            neighbours.push((neighbour_node.index(), weight));
            neighbour_string = format!("{} {}:{}", neighbour_string, neighbour_node.index(), weight);
        }
        let temp_string = format!("{} {}{}", node_index, base, neighbour_string);
        write_string = format!("{}\n{}", write_string, temp_string.clone());
    };
    write_string_to_newfile(&write_path, &write_string);
}

fn load_the_graph (file_name: String) -> Graph<u8, i32, Directed, usize> {
    let mut edge_capacity = 0;
    let mut max_node_index = 0;
    let mut node_loader: Vec<(usize, u8, Vec<(usize, i32)>)> = vec![];
    // check if available, populate node_edge_list from file
    if check_file_availability(&file_name, INTERMEDIATE_PATH) == true {
        // read the file
        let read_path = format!("{}/{}", INTERMEDIATE_PATH, file_name);
        // first popopulate the node list
        for line in read_to_string(&read_path).unwrap().lines() {
            let line_parts: Vec<&str> = line.split(" ").collect();
            // node definition
            if line_parts.len() >= 2 {
                let mut node_index = 0;
                let mut base = 0;
                let mut neighbours = vec![];
                let mut iter_index = 0;
                let mut line_part_iterator = line_parts.iter();
                while let Some(line_part) = line_part_iterator.next() {
                    if iter_index == 0 {
                        node_index = line_part.parse::<usize>().unwrap();
                        if node_index > max_node_index {
                            max_node_index = node_index;
                        }
                    }
                    else if iter_index == 1 {
                        base = line_part.parse::<u8>().unwrap();
                    }
                    else {
                        let neighbour_weight: Vec<&str> = line_part.split(":").collect();
                        let neighbour = neighbour_weight[0].parse::<usize>().unwrap();
                        let weight = neighbour_weight[1].parse::<i32>().unwrap();
                        neighbours.push((neighbour, weight));
                        edge_capacity += 1;
                    }
                    iter_index += 1;
                }
                node_loader.push((node_index, base, neighbours.clone()));
            }
        }
    }
    else {
        println!("file not available, skipping");
        return Graph::default();
    }
    // make the graph from the vector
    let mut graph: Graph<u8, i32, Directed, usize> = Graph::with_capacity(max_node_index + 1, edge_capacity);
    
    // add the nodes first
    for index in 0..max_node_index + 1 {
        match node_loader.iter().position(|r| r.0 == index) {
            Some(x) => {
                graph.add_node(node_loader[x].1);
            },
            None => {
                println!("{} not available", index);
                //graph.add_node(0);
            }
        }
    }
    // add the edges
    for node in node_loader{
        for edge in node.2 {
            graph.add_edge(NodeIndex::new(node.0), NodeIndex::new(edge.0), edge.1);
        }
    }
    // show graph 
    let write_string = format!("{}", Dot::new(&graph.map(|_, n| (*n) as char, |_, e| *e)));
    let write_path = format!("./result/loaded_graph.txt");
    write_string_to_newfile(&write_path, &write_string);
    graph
}

pub fn get_parallel_bases_from_file(file_path: &String, required_pos: usize) -> String {
    // open the file
    let f = File::open(&file_path).unwrap();
    let mut reader = BufReader::new(f);
    let mut buffer = String::new();
    let mut current_pos = 0;
    loop {
        buffer.clear();
        reader.read_line(&mut buffer).unwrap();
        if current_pos == required_pos {
            break;
        }
        current_pos += 1;
    }
    buffer
}

fn check_file_availability (file_name: &str, search_path: &str) -> bool {
    let temp_path_string = format!("{}/{}", search_path, file_name);
    let path = std::path::Path::new(&temp_path_string);
    let prefix = path.parent().unwrap();
    // check if dir is avaiable
    let file_available = match read_dir(prefix) {
        Ok(_) => {
            // check if file is available
            if path.exists() {
                true
            }
            else {
                false
            }
        },
        Err(_) => {false},
    };
    file_available
}

fn reverse_complement_filter_and_rearrange_subreads (original_subreads: &Vec<String>) -> Vec<String> {
    let mut seqvec: Vec<String> = vec![];
    //reverse complement every other line
    let mut index = 0;
    for seq in original_subreads {
        if index % 2 != 0 {
            let mut tempseq: Vec<char> = vec![];
            let iterator = seq.chars().rev().into_iter();
            for char in iterator{
                tempseq.push(match char {
                    'A' => 'T',
                    'C' => 'G',
                    'G' => 'C',
                    'T' => 'A',
                    _ => ' ',
                });
            }
            seqvec.push(tempseq.iter().cloned().collect::<String>());
        }
        else {
            seqvec.push((*seq.clone()).to_string());
        }
        index += 1;
    }
    //get rid of the last incomplete reading
    //seqvec.pop();
    //sort the vector by size
    seqvec.sort_by_key(|seq| seq.len());
    //drop the sequences which are > 1.8x median size
    let median_size: f32 = seqvec[(seqvec.len() / 2) - 1].len() as f32;
    let mut drop_index = seqvec.len();
    for index in (seqvec.len() / 2)..(seqvec.len() - 1) {
        if seqvec[index].len() as f32 > (median_size * 1.5) {
            drop_index = index;
            break;
        }
    }
    for _ in drop_index..seqvec.len() {
        seqvec.pop();
    }
    // rearrange the seq vector median first and rest according mediand size difference
    seqvec.sort_by(|a, b| ((a.len() as f32 - median_size).abs()).partial_cmp(&(b.len() as f32 - median_size).abs()).unwrap());
    seqvec
}

fn get_the_subreads_by_name_sam (full_name: &String) -> Vec<String> {
    let mut subread_vec: Vec<String> = vec![];
    let mut split_text_iter = (full_name.split("/")).into_iter();
    let file_name = split_text_iter.next().unwrap();
    let required_id = split_text_iter.next().unwrap().parse::<usize>().unwrap();
    let path = format!("{}{}{}", DATA_PATH.to_string(), file_name, ".subreads.sam".to_string());
    if file_name.eq(&"m64125_201017_124255".to_string()) {
        return subread_vec;
    }
    let file_position = read_index_file_for_sam (&file_name.to_string(), required_id);
    // file stuff init
    let mut count = 0;
    let f = File::open(&path).unwrap();
    let mut reader = BufReader::new(f);
    let mut buffer = String::new();
    reader.seek(SeekFrom::Start(file_position as u64)).expect("");
    loop {
        buffer.clear();
        reader.read_line(&mut buffer).unwrap();
        // split it to find the id
        let mut temp_split_iter = (buffer.split("/")).into_iter();
        temp_split_iter.next();
        let current_id;
        match temp_split_iter.next().unwrap().parse::<usize>() {
            Ok(x) => {current_id = x;},
            Err(_) => {break;},
        }
        if current_id != required_id {
            break;
        }
        else {
            // write code to extract the sequence and add to subread_vec
            let mut data_split_iter = (buffer.split("\t")).into_iter();
            println!("{}", data_split_iter.next().unwrap());
            for _ in 0..8 {data_split_iter.next();}
            subread_vec.push(data_split_iter.next().unwrap().to_string());
            count += 1;
        }
    }
    println!("count = {}", count);
    drop(reader);
    drop(buffer);
    subread_vec
}

pub fn read_index_file_for_sam (file_name: &String, read_name: usize) -> usize {
    // get the file location from index file if no index found make one
    println!("Reading index file {}{}", file_name, ".cui".to_string());
    let index_path = format!("result/{}.cui", file_name);
    let f;
    f = match File::open(&index_path) {
        Ok(x) => {x},
        Err(_) => {make_index_file_for_sam(&file_name.to_string())},
    };
    let mut reader = BufReader::new(f);
    let mut buffer = String::new();
    
    // get the file length
    reader.seek(SeekFrom::End(0)).expect("");
    let end_file_pos = reader.stream_position().unwrap();
    
    //go back to start
    reader.seek(SeekFrom::Start(0)).expect("");

    // go through and find index
    // jump to the middle of the file **binary searching**
    let mut section_start_file_pos = 0;
    let mut section_end_file_pos = end_file_pos;
    let mut current_file_pos = (section_start_file_pos + section_end_file_pos) / 2;
    let mut required_position = 0;
    loop {
        reader.seek(SeekFrom::Start(current_file_pos)).expect("");
        // get rid of the half line
        buffer.clear();
        reader.read_line(&mut buffer).unwrap();
        // the required line
        buffer.clear();
        reader.read_line(&mut buffer).unwrap();
        // split it to find the id
        let mut temp_split_iter = (buffer.split("\t")).into_iter();
        
        let current_id;
        match temp_split_iter.next().unwrap().parse::<usize>() {
            Ok(x) => {current_id = x;},
            Err(_) => {break;},
        }
        let required_position_string = temp_split_iter.next().unwrap().replace("\n", "");
        required_position = required_position_string.parse::<usize>().unwrap();
        //println!("curr: {}", current_id);
        // jumping 
        if read_name == current_id {
            break;
        }
        else if read_name > current_id {
            section_start_file_pos = current_file_pos;
        }
        else {
            section_end_file_pos = current_file_pos;
        }
        current_file_pos = (section_start_file_pos + section_end_file_pos) / 2;
    }
    required_position
}

pub fn make_index_file_for_sam (file_name: &String) -> File {
    println!("Making index file for {}{}", file_name, ".subreads.sam".to_string());
    let path = format!("{}{}{}", DATA_PATH.to_string(), file_name, ".subreads.sam".to_string());
    let write_path = format!("result/{}.cui", file_name);
    // get the file name and load it
    // file stuff init
    let f = File::open(&path).unwrap();
    let mut reader = BufReader::new(f);
    let mut buffer = String::new();

    // get the file length
    reader.seek(SeekFrom::End(0)).expect("");
    let end_file_pos = reader.stream_position().unwrap();
    
    //go back to start
    reader.seek(SeekFrom::Start(0)).expect("");
    
    // go through the file saving the sequence indices
    let mut write_string: String = "".to_string();
    let mut current_position;
    let mut current_ccs: usize;
    let mut prev_ccs: usize = 0;
    let mut index = 0;
    let mut break_count = 0;
    loop {
        // get the next line if available
        buffer.clear();
        match reader.read_line(&mut buffer) {
            Ok(_) => {if buffer.len() == 0 {break;}},
            Err(_) => {break;},
        }
        // split it to find the id
        let mut temp_split_iter = (buffer.split("/")).into_iter();
        temp_split_iter.next();
        match temp_split_iter.next() {
            Some(x) => {current_ccs = match x.parse::<usize>() { Ok(val) => {val}, Err(_) => {0}};},
            None => {continue;},
        }
        // add to the pos_read_name if different
        if current_ccs != prev_ccs {
            break_count = 0;
            prev_ccs = current_ccs;
            current_position = reader.stream_position().unwrap();
            write_string = format!("{}\n{}\t{}", write_string, current_ccs, current_position);
            index += 1;
            // display progress and write the current data
            if index % 100000 == 0 {
                println!("Progress {}%", (current_position * 100) / end_file_pos);
                write_string_to_file(&write_path, &write_string);
                write_string = "".to_string();
            }
        }
        else {
            break_count += 1;
            if break_count > 10000 {
                break;
            }
        }
    }
    // write the rest
    write_string_to_file(&write_path, &write_string);
    File::open(&write_path).unwrap()
}

fn get_the_subreads_by_name_bam (error_chr: &String, error_pos: usize, full_name: &String) -> Vec<String> {
    let subread_vec: Vec<String> = vec![];
    let mut split_text_iter = (full_name.split("/")).into_iter();
    let file_name = split_text_iter.next().unwrap();
    let required_id = split_text_iter.next().unwrap().parse::<i64>().unwrap();
    let path = format!("{}{}{}", DATA_PATH.to_string(), file_name, ".subreads.mapped.bam".to_string());
    if file_name.eq(&"m64125_201017_124255".to_string()) {
        return subread_vec;
    }
    let mut bam_reader = BamIndexedReader::from_path(path).unwrap();
    bam_reader.fetch((error_chr, error_pos as i64, error_pos as i64 + 1)).unwrap();
    let mut index = 0;
    for read in bam_reader.records() {
        let readunwrapped = read.unwrap();
        // get the required name
        let read_name = String::from_utf8(readunwrapped.qname().to_vec()).expect("");
        let mut split_read_name_iter = (read_name.split("/")).into_iter();
        split_read_name_iter.next().unwrap();
        let css_id = split_read_name_iter.next().unwrap().parse::<i64>().unwrap();
        if required_id != css_id {
            continue;
        }

        println!("readname {}", read_name);
        if readunwrapped.seq_len() < 5 {
            continue;
        }
        index += 1;
    }
    println!("new count = {}", index);
    
    subread_vec
}

fn check_the_scores_and_change_alignment (seqvec: Vec<String>, pacbio_consensus: &String) -> Vec<String> {
    let mut forward_score = 0;
    let mut backward_score = 0;
    // make the pacbio orientation files
    let pacbio_forward = pacbio_consensus.as_bytes().to_vec();
    let pacbio_backward;
    let mut tempseq: Vec<char> = vec![];
    let iterator = pacbio_consensus.chars().rev().into_iter();
    for char in iterator{
        tempseq.push(match char {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            _ => ' ',
        });
    }
    pacbio_backward = tempseq.iter().cloned().collect::<String>().as_bytes().to_vec();
    // check the forward scores for 2 sequences
    for seq in &seqvec {
        let score_func = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let k = 8; // kmer match length
        let w = 20; // Window size for creating the band
        let mut aligner = BandedDP::new(-5, -1, score_func, k, w);
        let alignment = aligner.local(&pacbio_forward, &seq.as_bytes().to_vec());
        let score = alignment.score;
        println!("forward score: {}", score);
        forward_score += score;
        break;
    }
    // check the backward scores for 2 sequences
    for seq in &seqvec {
        let score_func = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let k = 8; // kmer match length
        let w = 20; // Window size for creating the band
        let mut aligner = BandedDP::new(-5, -1, score_func, k, w);
        let alignment = aligner.local(&pacbio_backward, &seq.as_bytes().to_vec());
        let score = alignment.score;
        println!("backward score: {}", score);
        backward_score += score;
        break;
    }
    if forward_score < SKIP_SCORE && backward_score < SKIP_SCORE {
        return vec![];
    }
    else if backward_score > forward_score {
        println!("Scores are too low, inverting sequences.");
        let mut seqvec2 = vec![];
        //reverse complement every line
        for seq in &seqvec {
            let mut tempseq: Vec<char> = vec![];
            let iterator = seq.chars().rev().into_iter();
            for char in iterator{
                tempseq.push(match char {
                    'A' => 'T',
                    'C' => 'G',
                    'G' => 'C',
                    'T' => 'A',
                    _ => ' ',
                });
            }
            seqvec2.push(tempseq.iter().cloned().collect::<String>());
        }
        return seqvec2;
    }
    else {
        return seqvec;
    }
}

fn get_corrosponding_seq_name_location_quality_from_bam (error_pos: usize, error_chr: &String, base_change: &char) -> Vec<(String, String, u8, usize, Vec<f32>)> {
    let mut seq_name_qual_and_errorpos: Vec<(String, String, u8, usize, Vec<f32>)> = vec![]; // seq name qual errorpos 4 sn values
    let path = &READ_BAM_PATH;
    let mut bam_reader = BamIndexedReader::from_path(path).unwrap();
    bam_reader.fetch((error_chr, error_pos as i64, error_pos as i64 + 1)).unwrap();
    'read_loop: for read in bam_reader.records() {
        let readunwrapped = read.unwrap();
        // get the sn tag info
        for i in readunwrapped.aux_iter() {
            println!("{:?}", i.unwrap().0);
        };

        let mut sn_array= vec![0.0];
        readunwrapped.aux(b"sn").unwrap();
        // if let Ok(Aux::ArrayFloat(array)) = readunwrapped.aux(b"sn") {
        //     sn_array = array.iter().collect::<Vec<_>>();
        // }
        // else {
        //     panic!("Could not sn data");
        // }
        let mut read_index = 0;
        let read_name = String::from_utf8(readunwrapped.qname().to_vec()).expect("");
        let read_vec = readunwrapped.seq().as_bytes().to_vec();
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
                    if (current_ref_pos + temp_int > error_pos)
                        && (current_ref_pos <= error_pos + 1) {
                        (read_index, _) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, error_pos, 1);
                        if &(read_vec[read_index] as char) == base_change && ('X' != *base_change) {
                            break;
                        }
                        else if 'X' == *base_change {
                            break;
                        }
                        else {
                            continue 'read_loop;
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
                    if (current_ref_pos + temp_int >= error_pos)
                        && (current_ref_pos <= error_pos + 1) {
                        //(_, _) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, error_pos, 1);
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
                        //let (_, _) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, error_pos, 1);
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
        seq_name_qual_and_errorpos.push((read_string.clone(), read_name.clone(), readunwrapped.qual()[read_index], read_index, sn_array));
    }
    drop(bam_reader);
    seq_name_qual_and_errorpos
}

pub fn write_string_to_file (file_name: &String, input_string: &String) {
    let path = std::path::Path::new(&file_name);
    let prefix = path.parent().unwrap();
    create_dir_all(prefix).unwrap();
    let mut file = OpenOptions::new().create(true).append(true).open(file_name).unwrap();
    write!(file, "{}", input_string).expect("result file cannot be written");
}

pub fn get_zoomed_graph_section (normal_graph: &Graph<u8, i32, Directed, usize>, focus_node: &usize)-> String {
    let mut graph_section= "".to_string();
    let normal_dot = format!("{:?}", Dot::new(&normal_graph.map(|_, n| (*n) as char, |_, e| *e)));
    let displaying_nodes: Vec<usize> = find_neighbouring_indices (NUM_OF_ITER_FOR_ZOOMED_GRAPHS, *focus_node, normal_graph);
    let mut graph_section_nodes: String = "".to_string();
    let mut graph_section_edges: String = "".to_string();
    //find the position in the dot file and add to the graph section
    for node in &displaying_nodes {
        //get the nodes from the dot file
        match normal_dot.find(&format!(" {} [", node)) {
            Some(start) => {
                let mut char_seq = vec![];
                let mut end = start;
                loop {
                    char_seq.push(normal_dot.chars().nth(end).unwrap());
                    if normal_dot.chars().nth(end).unwrap() == '\n' {
                        break;
                    }
                    end += 1;
                }
                //add from start to end to the graph_section string
                graph_section_nodes = char_seq.iter().collect::<String>();
            },
            None => {}
        }
        graph_section = format!("{}{}", graph_section, graph_section_nodes).to_string();
        graph_section_nodes = "".to_string();
    }
    for node in &displaying_nodes {
        //get the edges from the dot file
        let edge_entries: Vec<usize> = normal_dot.match_indices(&format!(" {} ->", node)).map(|(i, _)|i).collect();
        for edge in edge_entries {
            let mut char_seq = vec![];
                let mut end = edge;
                loop {
                    char_seq.push(normal_dot.chars().nth(end).unwrap());
                    if normal_dot.chars().nth(end).unwrap() == '\n' {
                        break;
                    }
                    end += 1;
                }
                //add from start to end to the graph_section string
                graph_section_edges = format!("{}{}", graph_section_edges, char_seq.iter().collect::<String>()).to_string();
        }
        graph_section = format!("{}{}", graph_section, graph_section_edges).to_string();
        graph_section_edges = "".to_string();
    }
    //modifying the section graph with the highlight
    graph_section = modify_dot_graph_with_highlight (graph_section, focus_node);
    //make it a dot graph
    graph_section = format!("digraph {{\n{} }}", graph_section);
    graph_section
}


fn find_neighbouring_indices (num_of_iterations: usize, focus_node: usize, graph: &Graph<u8, i32, Directed, usize> ) -> Vec<usize> {
    let mut indices: Vec<usize> = vec![];
    if num_of_iterations <= 0 {
        return indices;
    }
    let mut immediate_neighbours = graph.neighbors_undirected(NodeIndex::new(focus_node));
    while let Some(neighbour_node) = immediate_neighbours.next() {
        if !indices.contains(&neighbour_node.index()){
            indices.push(neighbour_node.index());
        }
        let obtained_indices = find_neighbouring_indices(num_of_iterations - 1, neighbour_node.index(), graph);
        for obtained_index in obtained_indices {
            if !indices.contains(&obtained_index){
                indices.push(obtained_index);
            }
        }
    }
    indices 
}

fn modify_dot_graph_with_highlight (mut dot: String, focus_node: &usize) -> String {
    match dot.find(&format!(" {} [", focus_node)) {
        Some(mut x) => {
            while dot.chars().nth(x).unwrap() != '[' {
                x += 1;
            }
            dot.replace_range(x + 12..x + 12,"Target node");
        },
        None => {}
    };
    dot
}

fn read_fai_get_ref_context (start_pos: usize, length: usize, chromosone: &String, reader: &mut faidx::Reader) -> String {
    reader.fetch_seq_string(chromosone, start_pos, start_pos + length - 1).unwrap()
}

pub fn get_required_start_end_positions_from_read (section_length: usize, current_ref_pos: usize, current_read_pos: usize, required_pos: usize, required_len: usize) -> (usize, usize) {
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
