mod cli;
mod types;
mod io;
mod calculation;

use clap::Parser;
use crate::cli::{Args, Formula, InputMode};
use crate::types::{Calculation, FilterCriteria, SiteData};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

fn main() {
    let args = Args::parse();

    let (input_mode, files) = args.get_input_files().unwrap_or_else(|e| {
        eprintln!("{}", e);
        std::process::exit(1);
    });

    let criteria = FilterCriteria {
        min_depth: args.min_depth,
        min_freq: args.min_af,
        min_alt_reads: args.min_alt_reads,
        min_alt_rev_reads: args.min_alt_rev_reads,
    };

    // Inform user about applied filters
    eprintln!("\n--- Applying Filters ---");
    eprintln!("> Minimum Depth (DP): {}", criteria.min_depth);
    eprintln!("> Minimum Allele Freq (AF): {}", criteria.min_freq);
    eprintln!("> Minimum Alternate Reads (AD): {}", criteria.min_alt_reads);
    eprintln!("> Minimum Alt. Reverse Reads (ADR): {}", criteria.min_alt_rev_reads);
    eprintln!("------------------------\n");

    let formula = match args.formula {
        Formula::Fst => Calculation::Fst,
        Formula::Gst => Calculation::Gst,
        Formula::Nei => Calculation::Nei,
        Formula::Chord => Calculation::Chord,
        Formula::BrayCurtis => Calculation::BrayCurtis,
        Formula::JostD => Calculation::JostD,
        Formula::Reynolds => Calculation::Reynolds,
        Formula::Rogers => Calculation::Rogers,
    };

    args.configure_thread_pool();

    // Read reference if provided
    let reference_seq = args.reference.as_ref().map(|path| {
        eprintln!("Reading reference FASTA...");
        io::fasta::read_reference_sequence(path.to_str().unwrap_or_default()).unwrap_or_else(|e| {
            eprintln!("Error reading reference file: {}", e);
            std::process::exit(1);
        })
    });

    // Read input files
    eprintln!("Reading and processing input files...");
    let (all_positions, mut variants_by_sample) = match input_mode {
        InputMode::Vcf => io::vcf::read_vcf_files(&files, &criteria),
        InputMode::Table => io::csv::read_csv_files(&files, &criteria),
    };

    // FASTA fallback for reference allele (only when reference is provided)
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

    eprintln!(
        "Found {} samples and {} polymorphic sites after filtering.",
        variants_by_sample.len(),
        all_positions.len()
    );

    // Prepare for calculation
    let samples: Vec<String> = {
        let mut s: Vec<String> = variants_by_sample.keys().cloned().collect();
        s.sort();
        s
    };
    let sample_pairs = calculation::generate_sample_pairs(&samples);

    // Pairwise calculations
    let num_loci = all_positions.len();
    let output_path_str = args.output.to_str().unwrap_or("output");

    eprintln!(
        "Computing {} distances for {} pairs...",
        format!("{:?}", args.formula).to_lowercase(),
        sample_pairs.len()
    );

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
            let normalize = args.normalize;

            let dist = match formula {
                Calculation::Fst => calculation::fst::calculate_fst_for_pair(data1, data2, &all_positions, normalize, num_loci),
                Calculation::Gst => calculation::gst::calculate_gst_for_pair(data1, data2, &all_positions),
                Calculation::Nei => calculation::nei::calculate_nei_distance_for_pair(data1, data2, &all_positions),
                Calculation::Chord => calculation::chord::calculate_chord_distance_for_pair(data1, data2, &all_positions, normalize, num_loci),
                Calculation::BrayCurtis => calculation::bray_curtis::calculate_bray_curtis_for_pair(data1, data2, &all_positions, normalize, num_loci),
                Calculation::JostD => calculation::jost_d::calculate_jost_d_for_pair(data1, data2, &all_positions, normalize, num_loci),
                Calculation::Reynolds => calculation::reynolds::calculate_reynolds_distance_for_pair(data1, data2, &all_positions),
                Calculation::Rogers => calculation::rogers::calculate_rogers_distance_for_pair(data1, data2, &all_positions, num_loci),
            };
            pb.inc(1);
            ((i, j), dist)
        })
        .collect();

    pb.finish_with_message("Calculation complete.");

    // Write output
    let distance_matrix = calculation::create_distance_matrix(&results, samples.len());

    // Auto-detect output delimiter from extension
    let use_tab = output_path_str.ends_with(".tsv") || output_path_str.ends_with(".tab");

    match io::csv::write_distance_matrix(output_path_str, &distance_matrix, &samples, use_tab) {
        Ok(_) => eprintln!("Distance matrix written to {}", output_path_str),
        Err(e) => {
            eprintln!("Error writing output file: {}", e);
            std::process::exit(1);
        }
    }
}
