mod cli;
mod types;
mod io;
mod calculation;

use crate::types::{Calculation, FilterCriteria, SiteData};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

fn main() {
    let args = cli::Cli::parse();

    let (input_mode, files) = args.get_input_files().unwrap_or_else(|e| {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    });

    // Require reference for table mode
    if matches!(input_mode, cli::InputMode::Table) && args.reference.is_none() {
        eprintln!("Error: --reference is required for table mode.");
        std::process::exit(1);
    }

    let criteria = FilterCriteria {
        min_depth: args.min_depth,
        min_freq: args.min_af,
        min_alt_reads: args.min_alt_reads,
        min_alt_rev_reads: args.min_alt_rev_reads,
    };

    println!("\n--- Applying Filters ---");
    println!("> Minimum Depth (DP): {}", criteria.min_depth);
    println!("> Minimum Allele Freq (AF): {}", criteria.min_freq);
    println!("> Minimum Alternate Reads (AD): {}", criteria.min_alt_reads);
    println!("> Minimum Alt. Reverse Reads (ADR): {}", criteria.min_alt_rev_reads);
    println!("------------------------\n");

    let formula = match args.formula.as_str() {
        "fst" => Calculation::Fst,
        "gst" => Calculation::Gst,
        "nei" => Calculation::Nei,
        "chord" => Calculation::Chord,
        "bray-curtis" => Calculation::BrayCurtis,
        "jost_d" => Calculation::JostD,
        "reynolds" => Calculation::Reynolds,
        "rogers" => Calculation::Rogers,
        _ => unreachable!(),
    };

    args.configure_thread_pool();

    // Read reference (if provided)
    let reference_seq = if let Some(ref ref_path) = args.reference {
        println!("Reading reference FASTA...");
        let path_str = ref_path.to_str().unwrap_or("");
        Some(io::fasta::read_reference_sequence(path_str).unwrap_or_else(|e| {
            eprintln!("Error reading reference file {}: {}", path_str, e);
            std::process::exit(1);
        }))
    } else {
        None
    };

    println!("Reading and processing input files...");
    let (all_positions, mut variants_by_sample) = match input_mode {
        cli::InputMode::Vcf => io::vcf::read_vcf_files(&files, &criteria),
        cli::InputMode::Table => io::csv::read_csv_files(&files, &criteria),
    };

    // FASTA fallback for reference allele if needed
    if let Some(ref ref_seq) = reference_seq {
        for pos in &all_positions {
            for sample_data in variants_by_sample.values_mut() {
                let site = sample_data.entry(*pos).or_insert_with(SiteData::default);
                if site.reference_allele.is_empty() {
                    if let Some(ref_char) = ref_seq.get(pos - 1) {
                        site.reference_allele = ref_char.to_string();
                    }
                }
            }
        }
    }

    println!("Found {} samples and {} polymorphic sites after filtering.", variants_by_sample.len(), all_positions.len());

    let samples: Vec<String> = {
        let mut s: Vec<String> = variants_by_sample.keys().cloned().collect();
        s.sort();
        s
    };
    let sample_pairs = calculation::generate_sample_pairs(&samples);
    let num_loci = all_positions.len();

    let pb = ProgressBar::new(sample_pairs.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({percent}%)")
            .unwrap()
            .progress_chars("#>-"),
    );

    let results: Vec<((usize, usize), f64)> = sample_pairs
        .par_iter()
        .map(|&(i, j)| {
            let data1 = variants_by_sample.get(&samples[i]).unwrap();
            let data2 = variants_by_sample.get(&samples[j]).unwrap();

            let dist = match formula {
                Calculation::Fst => calculation::fst::calculate_fst_for_pair(data1, data2, &all_positions, args.normalize, num_loci),
                Calculation::Gst => calculation::gst::calculate_gst_for_pair(data1, data2, &all_positions),
                Calculation::Nei => calculation::nei::calculate_nei_distance_for_pair(data1, data2, &all_positions),
                Calculation::Chord => calculation::chord::calculate_chord_distance_for_pair(data1, data2, &all_positions, args.normalize, num_loci),
                Calculation::BrayCurtis => calculation::bray_curtis::calculate_bray_curtis_for_pair(data1, data2, &all_positions, args.normalize, num_loci),
                Calculation::JostD => calculation::jost_d::calculate_jost_d_for_pair(data1, data2, &all_positions, args.normalize, num_loci),
                Calculation::Reynolds => calculation::reynolds::calculate_reynolds_distance_for_pair(data1, data2, &all_positions),
                Calculation::Rogers => calculation::rogers::calculate_rogers_distance_for_pair(data1, data2, &all_positions, num_loci),
            };
            pb.inc(1);
            ((i, j), dist)
        })
        .collect();

    pb.finish_with_message("Calculation complete.");

    let distance_matrix = calculation::create_distance_matrix(&results, samples.len());
    let output_str = args.output.to_str().unwrap_or("output");
    match io::csv::write_distance_matrix(output_str, &distance_matrix, &samples) {
        Ok(_) => println!("Distance matrix successfully written to {}", output_str),
        Err(e) => eprintln!("Error writing output file: {}", e),
    }
}
