#![allow(dead_code)]
mod alignment;
mod generator;
mod misc;
mod quality;

use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use std::io::SeekFrom;
use pprof;

use petgraph::Directed;
use petgraph::Graph;
use quality::topology_cut::base_quality_score_calculation;
use quality::topology_cut::get_parallel_nodes_with_topology_cut;
use alignment::poarustbio::Aligner;
use alignment::pairwise::pairwise;
use misc::get_error_bases_from_himut_vcf;
use misc::get_required_start_end_positions_from_read;
use rust_htslib::bam::{Read as BamRead, IndexedReader as BamIndexedReader};
use misc::get_zoomed_graph_section;
use misc::write_string_to_file;

const SEED: u64 = 2;
const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = 0;
const MATCH: i32 = 2;
const MISMATCH: i32 = -2;
const DATA_PATH: &str = "/data1/hifi_consensus/try2/";
const READ_BAM_PATH: &str = "/data1/hifi_consensus/try2/merged.bam";
fn main() {
    let guard = pprof::ProfilerGuardBuilder::default().frequency(1000).blocklist(&["libc", "libgcc", "pthread", "vdso"]).build().unwrap();
    pipeline_redo_poa_get_topological_quality_score();
    if let Ok(report) = guard.report().build() {
        let file = File::create("result/flamegraph.svg").unwrap();
        report.flamegraph(file).unwrap();
    };
    drop(guard);
}

fn pipeline_redo_poa_get_topological_quality_score () {
    // get the error locations
    let error_locations = get_error_bases_from_himut_vcf (); //chromosone, location, ref allele, alt allele
    // go through the error locations
    for error_location in error_locations {
        // start after this error location, 
        let skip_location = 13405643;
        let skip_chromosone = "chr1";
        if (error_location.0 == skip_chromosone) && (error_location.1 < skip_location) {
            continue;
        }
        println!("Error position {}:{} ref allele: {} alt allele: {}", error_location.0, error_location.1, error_location.2, error_location.3);
        // find the ccs which are in that error
        let seq_name_qual_and_errorpos_vec = get_corrosponding_seq_name_location_quality_from_bam(error_location.1, &error_location.0, &error_location.3);
        for seq_name_qual_and_errorpos in seq_name_qual_and_errorpos_vec {
            println!("Processing ccs file: {}", seq_name_qual_and_errorpos.1);
            // find the subreads of that ccs
            let mut sub_reads = get_the_subreads_by_name(&seq_name_qual_and_errorpos.1);
            // skip if no subreads, errors and stuff
            if sub_reads.len() == 0 {
                continue;
            }
            // reverse the sub reads if score is low
            sub_reads = check_the_scores_and_change_alignment(sub_reads, &seq_name_qual_and_errorpos.0);
            
            sub_reads.insert(0, seq_name_qual_and_errorpos.0.clone());
            // do poa with the read and subreads, get the poa and consensus
            let mut sequence_number: usize = 0;
            let mut aligner = Aligner::new(MATCH, MISMATCH, GAP_OPEN, &sub_reads[0].as_bytes().to_vec());
            for sub_read in &sub_reads {
                if sequence_number != 0 {
                    aligner.global(&sub_read.as_bytes().to_vec()).add_to_graph();
                }
                sequence_number += 1;
                println!("Sequence {} processed", sequence_number);
            }
            let (calculated_consensus, calculated_topology) = aligner.poa.consensus(); //just poa
            let calculated_graph: &Graph<u8, i32, Directed, usize> = aligner.graph();
            // check if fixed (check the consensus location which is matched to the read error location)
            let position = get_redone_consensus_error_position(seq_name_qual_and_errorpos.0, &calculated_consensus, seq_name_qual_and_errorpos.3);
            // calculate the quality score of the location
            let skip_nodes: Vec<usize> = calculated_topology[0 .. position + 1].to_vec();
            let target_node_parent = Some(calculated_topology[position - 1]);
            let target_node_child = Some(calculated_topology[position + 1]);
            
            let (parallel_nodes, parallel_num_incoming_seq, _) = get_parallel_nodes_with_topology_cut (skip_nodes, sequence_number,  calculated_topology[position], target_node_parent, target_node_child, calculated_graph);
            let (calculated_quality_score, _, parallel_bases, _) = base_quality_score_calculation (sequence_number, parallel_nodes, parallel_num_incoming_seq, calculated_consensus[position], calculated_graph);
            let write_string = format!("Error position {}:{} ref allele: {} alt allele: {}\nPacbio base: \t{} quality: {}\nCalculated base: \t{} quality: {}\nParallel Bases: ACGT:{:?}\n\n", error_location.0, error_location.1, error_location.2, error_location.3, error_location.3, seq_name_qual_and_errorpos.2, calculated_consensus[position] as char, calculated_quality_score, parallel_bases);
            write_string_to_file("result/quality.txt", &write_string);
            let write_string = format!("{}\n{}\n\n", write_string, get_zoomed_graph_section(calculated_graph, &calculated_topology[position]));
            write_string_to_file("result/graph.txt", &write_string);
            break;
        }
        break;
    }
}

fn get_redone_consensus_error_position (pacbio_consensus: String, calculated_consensus: &Vec<u8>, pacbio_error_position: usize) -> usize {
    let mut consensus_match_invalid_indices: Vec<usize> = vec![];
    let pacbio_consensus_vec: Vec<u8> = pacbio_consensus.bytes().collect();
    let mut aligned_pacbio_scores_vec: Vec<usize> = vec![];
    let mut aligned_pacbio_bases:Vec<u8> = vec![];
    let (alignment, _) = pairwise(&calculated_consensus, &pacbio_consensus_vec, MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND);
    let mut pacbio_index = 0;
    let mut calc_index = 0;
    let mut calc_error_position: usize = 0;
    for op in alignment {
        match op as char {
            'm' => {
                pacbio_index += 1;
                calc_index += 1;
            },
            's' => {
                aligned_pacbio_bases.push(pacbio_consensus_vec[pacbio_index]);
                consensus_match_invalid_indices.push(calc_index);
                pacbio_index += 1;
                calc_index += 1;
            },
            'd' => {
                pacbio_index += 1;
            },
            'i' => {
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

fn check_the_scores_and_change_alignment (seqvec: Vec<String>, pacbio_consensus: &String) -> Vec<String> {
    let mut invert: bool = false;
    let mut seqvec2: Vec<String> = vec![];
    // check the scores for 3 sequences
    let mut index = 0;
    for seq in &seqvec {
        let (_, score) = pairwise(&pacbio_consensus.as_bytes().to_vec(), &seq.as_bytes().to_vec(), MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND);
        if score < 1000 {
            invert = true;
            break;
        }
        else if index > 3 {
            break;
        }
        index += 1;
    }
    if invert {
        println!("Scores are too low, inverting sequences.");
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
    }
    else {
        seqvec2 = seqvec;
    }
    seqvec2
}

fn get_the_subreads_by_name (full_name: &String) -> Vec<String> {
    let mut subread_vec: Vec<String> = vec![];
    let mut split_text_iter = (full_name.split("/")).into_iter();
    let file_name = split_text_iter.next().unwrap();
    let required_id = split_text_iter.next().unwrap().parse::<i64>().unwrap();
    let path = format!("{}{}{}", DATA_PATH.to_string(), file_name, ".subreads.sam".to_string());
    if file_name.eq(&"m64125_201017_124255".to_string()) {
        return subread_vec;
    }
    // file stuff init
    let f = File::open(&path).unwrap();
    let mut reader = BufReader::new(f);
    let mut buffer = String::new();

    // get the file length
    reader.seek(SeekFrom::End(0)).expect("");
    let end_file_pos = reader.stream_position().unwrap();

    // jump to the middle of the file
    let mut section_start_file_pos = 0;
    let mut section_end_file_pos = end_file_pos;
    let mut current_file_pos = (section_start_file_pos + section_end_file_pos) / 2;
    loop {
        reader.seek(SeekFrom::Start(current_file_pos)).expect("");
        // get rid of the half line
        buffer.clear();
        reader.read_line(&mut buffer).unwrap();
        // the required line
        buffer.clear();
        reader.read_line(&mut buffer).unwrap();
        // split it to find the id
        let mut temp_split_iter = (buffer.split("/")).into_iter();
        temp_split_iter.next();
        let current_id;
        match temp_split_iter.next().unwrap().parse::<i64>() {
            Ok(x) => {current_id = x;},
            Err(_) => {break;},
        }
        //println!("curr: {}", current_id);
        // jumping 
        if required_id == current_id {
            // when the id is found go back until you reach a different id
            let mut index = 1;
            loop {
                if current_file_pos < index * 10000 {
                    reader.seek(SeekFrom::Start(0)).expect("");
                    break;
                }
                reader.seek(SeekFrom::Start(current_file_pos - index * 100000)).expect("");
                // get rid of the half line
                buffer.clear();
                reader.read_line(&mut buffer).unwrap();
                // the required line
                buffer.clear();
                reader.read_line(&mut buffer).unwrap();
                // split it to find the id
                let mut temp_split_iter = (buffer.split("/")).into_iter();
                temp_split_iter.next();
                let current_id;
                match temp_split_iter.next().unwrap().parse::<i64>() {
                    Ok(x) => {current_id = x;},
                    Err(_) => {break;},
                }
                if current_id != required_id {
                    break;
                } 
                index += 1;
            }
            break;
        }
        else if required_id > current_id {
            section_start_file_pos = current_file_pos;
        }
        else {
            section_end_file_pos = current_file_pos;
        }
        current_file_pos = (section_start_file_pos + section_end_file_pos) / 2;
    }
    // count the number of entries and save the sequences
    let mut count = 0;
    let mut continue_count = 0;
    let continue_threshold = 3;
    loop {
        buffer.clear();
        reader.read_line(&mut buffer).unwrap();
        // split it to find the id
        let mut temp_split_iter = (buffer.split("/")).into_iter();
        temp_split_iter.next();
        let current_id;
        match temp_split_iter.next().unwrap().parse::<i64>() {
            Ok(x) => {current_id = x;},
            Err(_) => {break;},
        }
        if current_id != required_id {
            continue_count += 1;
            if continue_count > continue_threshold {
                break;
            }
        }
        else {
            // write code to extract the sequence and add to subread_vec
            let mut data_split_iter = (buffer.split("\t")).into_iter();
            for _ in 0..9 {data_split_iter.next();}
            subread_vec.push(data_split_iter.next().unwrap().to_string());
            count += 1;
        }
    }
    println!("count = {}", count);
    subread_vec
}

fn get_corrosponding_seq_name_location_quality_from_bam (error_pos: usize, error_chr: &String, base_change: &char) -> Vec<(String, String, u8, usize)> {
    let mut seq_name_qual_and_errorpos: Vec<(String, String, u8, usize)> = vec![];
    let path = &READ_BAM_PATH;
    let mut bam_reader = BamIndexedReader::from_path(path).unwrap();
    bam_reader.fetch((error_chr, error_pos as i64, error_pos as i64 + 1)).unwrap();
    'read_loop: for read in bam_reader.records() {
        let readunwrapped = read.unwrap();
        // get the data
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
                    if (current_ref_pos + temp_int >= error_pos)
                        && (current_ref_pos <= error_pos + 1) {
                        (read_index, _) = get_required_start_end_positions_from_read (temp_int, current_ref_pos, current_read_pos, error_pos, 1);
                        if &(read_vec[read_index] as char) == base_change {
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
        seq_name_qual_and_errorpos.push((read_string.clone(), read_name.clone(), readunwrapped.qual()[read_index], read_index));
    }
    seq_name_qual_and_errorpos
}
