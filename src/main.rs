mod cli;
mod types;
mod io;
mod calculation;

use clap::Parser;
use crate::cli::{Args, Formula, InputMode};
use crate::types::{Calculation, FilterCriteria, SiteData};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::process::ExitCode;

fn run() -> Result<(), String> {
    let args = Args::parse();

    let (input_mode, files) = args.get_input_files()?;

    let criteria = FilterCriteria {
        min_depth: args.min_depth,
        min_freq: args.min_af,
        min_alt_reads: args.min_alt_reads,
        min_alt_rev_reads: args.min_alt_rev_reads,
    };

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
    let reference = args.reference.as_ref().map(|path| {
        let path_str = path.to_str().unwrap_or_default();
        eprintln!("Reading reference FASTA...");
        io::fasta::read_reference(path_str)
    }).transpose().map_err(|e| format!("Reference error: {}", e))?;

    // Read input files
    eprintln!("Reading and processing input files...");
    let (all_positions, mut variants_by_sample) = match input_mode {
        InputMode::Vcf => io::vcf::read_vcf_files(&files, &criteria),
        InputMode::Table => io::csv::read_csv_files(&files, &criteria),
    };

    // FASTA fallback for reference allele (only when reference is provided)
    if let Some(ref genome) = reference {
        for pos in &all_positions {
            for sample_data in variants_by_sample.values_mut() {
                let site = sample_data.entry(*pos).or_insert_with(SiteData::default);
                if site.reference_allele.is_empty() {
                    if let Some(base) = io::fasta::get_base_at(genome, *pos) {
                        site.reference_allele = String::from(base as char);
                    }
                }
            }
        }
    }

    if variants_by_sample.is_empty() {
        return Err("No samples found after filtering. Check input files and filter thresholds.".to_string());
    }

    if all_positions.is_empty() {
        return Err("No polymorphic sites found after filtering.".to_string());
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

    if sample_pairs.is_empty() {
        return Err("Need at least 2 samples for pairwise distance calculation.".to_string());
    }

    let num_loci = all_positions.len();
    let output_path_str = args.output.to_str().ok_or("Invalid output path")?;

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
    let use_tab = output_path_str.ends_with(".tsv") || output_path_str.ends_with(".tab");

    io::csv::write_distance_matrix(output_path_str, &distance_matrix, &samples, use_tab)
        .map_err(|e| format!("Error writing output: {}", e))?;

    eprintln!("Distance matrix written to {}", output_path_str);
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("Error: {}", e);
            ExitCode::FAILURE
        }
    }
}
