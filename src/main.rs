

#![allow(dead_code)]
mod alignment;
mod generator;
mod misc;
mod quality;
use std::thread;
use crate::misc::calculate_deep_quality;
use rust_htslib::bcf::{self, Read as BcfRead};

const SEED: u64 = 2;
const GAP_OPEN: i32 = -2;
const GAP_EXTEND: i32 = 0;
const MATCH: i32 = 2;
const MISMATCH: i32 = -2;
const NTHREADS: usize = 64;
//m64125_201110_063134/0
fn main() {
    test_function();
    // make a vector to hold the children which are spawned.
    let mut children = vec![];
    let chromosone = "chr2";
    let total_start = 5_000_000;
    let total_end = 240_000_000;
    let one_thread_allocation = (total_end - total_start) / NTHREADS;
    //concancate_files();
    for i in 0..NTHREADS {
        // calculate my start and end locations
        // spin up another thread
        children.push(thread::spawn(move || {
            // calculate my start and end locations
            let start = total_start + one_thread_allocation * i;
            let end = start + one_thread_allocation;
            println!("Thread number {} started, {} from {} to {}..", chromosone, i, start, end);
            //calculate_deep_quality(chromosone, start, end, i);
        }));
    }
    for child in children {
        // wait for the thread to finish. Returns a result.
        let _ = child.join();
    }
}
 //to test vcf reading
fn test_function () {
    let mut vcf_reader = bcf::IndexedReader::from_path("/data1/phasstphase_test/potato/hifi/potato_deep.g.vcf.gz").expect("could not load index for vcf... looking for .csi file");
    match vcf_reader.fetch(0, 0 as u64, Some(10000)) {
        Ok(_) => {
            let mut total = 0;
            let mut hets = 0;
            for (i, _rec) in vcf_reader.records().enumerate() {
                total += 1;
                let rec = _rec.expect("cant unwrap vcf record");
                let pos = rec.pos();
                let alleles = rec.alleles();
                if alleles.len() > 2 {
                    continue; // ignore multi allelic sites
                }
                let reference = std::str::from_utf8(alleles[0]).expect("this really shouldnt fail");
                let alternative = std::str::from_utf8(alleles[1]). expect("this really shouldnt fail2");
                let rec_chrom = rec.rid().expect("could not unwrap vcf record id");
                let genotypes = rec.genotypes().expect("cant get genotypes");
                let genotype = genotypes.get(0); // assume only 1 and get the first one
                println!("{} {} {}", reference, alternative, genotypes.get(0));
            }
        },
        Err(_e) => (),
    }; // skip to chromosome for this thread
}

