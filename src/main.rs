mod cli;
mod types;
mod io;
mod calculation;

use crate::types::{Calculation, FilterCriteria, SiteData};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

fn main() {
    // 1. Parse command-line arguments
    let matches = cli::build_cli().get_matches();
    let (input_mode, files) = cli::get_input_files(&matches).unwrap_or_else(|e| {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    });

    let reference_path = matches.value_of("reference");
    let output_path = matches.value_of("output").unwrap();
    let normalize = matches.is_present("normalize");

    // Require reference for table mode
    if matches!(input_mode, cli::InputMode::Table) && reference_path.is_none() {
        eprintln!("Error: --reference is required for table mode.");
        std::process::exit(1);
    }

    // Parse filter arguments
    let criteria = FilterCriteria {
        min_depth: matches.value_of("min_depth").unwrap().parse().expect("min-depth must be an integer"),
        min_freq: matches.value_of("min_af").unwrap().parse().expect("min-af must be a float"),
        min_alt_reads: matches.value_of("min_alt_reads").unwrap().parse().expect("min-alt-reads must be an integer"),
        min_alt_rev_reads: matches.value_of("min_alt_rev_reads").unwrap().parse().expect("min-alt-rev-reads must be an integer"),
    };

    // Inform user about applied filters
    println!("\n--- Applying Filters ---");
    println!("> Minimum Depth (DP): {}", criteria.min_depth);
    println!("> Minimum Allele Freq (AF): {}", criteria.min_freq);
    println!("> Minimum Alternate Reads (AD): {}", criteria.min_alt_reads);
    println!("> Minimum Alt. Reverse Reads (ADR): {}", criteria.min_alt_rev_reads);
    println!("------------------------\n");

    let formula = match matches.value_of("formula").unwrap() {
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
    cli::configure_thread_pool(&matches);

    // 2. Read reference (if provided) and input files
    let reference_seq = if let Some(ref_path) = reference_path {
        println!("Reading reference FASTA...");
        Some(io::fasta::read_reference_sequence(ref_path).unwrap_or_else(|e| {
            eprintln!("Error reading reference file {}: {}", ref_path, e);
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

    // 3. FASTA Fallback for reference allele if needed
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

    // 4. Prepare for calculation
    let samples: Vec<String> = {
        let mut s: Vec<String> = variants_by_sample.keys().cloned().collect();
        s.sort();
        s
    };
    let sample_pairs = calculation::generate_sample_pairs(&samples);

    // 5. Perform pairwise calculations
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
            let sample1_name = &samples[i];
            let sample2_name = &samples[j];
            let data1 = variants_by_sample.get(sample1_name).unwrap();
            let data2 = variants_by_sample.get(sample2_name).unwrap();

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

    // 6. Write output matrix
    let distance_matrix = calculation::create_distance_matrix(&results, samples.len());
    match io::csv::write_distance_matrix(output_path, &distance_matrix, &samples) {
        Ok(_) => println!("Distance matrix successfully written to {}", output_path),
        Err(e) => eprintln!("Error writing output file: {}", e),
    }
}
