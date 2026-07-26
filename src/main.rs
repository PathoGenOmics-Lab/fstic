mod cli;
mod types;
mod io;
mod calculation;

use clap::Parser;
use crate::cli::{Args, Formula, InputMode};
use crate::types::FilterCriteria;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::process::ExitCode;
use std::time::Instant;

fn run() -> Result<(), String> {
    let args = Args::parse();
    let start = Instant::now();

    let (input_mode, files) = args.get_input_files()?;

    for f in &files {
        if !f.exists() {
            return Err(format!("Input file not found: {}", f.display()));
        }
    }

    let criteria = FilterCriteria {
        min_depth: args.min_depth,
        min_freq: args.min_af,
        min_alt_reads: args.min_alt_reads,
        min_alt_rev_reads: args.min_alt_rev_reads,
        pass_only: args.pass_only,
    };

    eprintln!("\n--- Applying Filters ---");
    eprintln!("> Minimum Depth (DP): {}", criteria.min_depth);
    eprintln!("> Minimum Allele Freq (AF): {}", criteria.min_freq);
    eprintln!("> Minimum Alternate Reads (AD): {}", criteria.min_alt_reads);
    eprintln!("> Minimum Alt. Reverse Reads (ADR): {}", criteria.min_alt_rev_reads);
    eprintln!("------------------------\n");

    args.configure_thread_pool();

    // Read reference if provided
    let reference = args.reference.as_ref().map(|path| {
        let path_str = path.to_str().unwrap_or_default();
        eprintln!("Reading reference FASTA...");
        io::fasta::read_reference(path_str)
    }).transpose().map_err(|e| format!("Reference error: {}", e))?;

    // Read input files
    eprintln!("Reading {} input files...", files.len());
    let (position_set, mut variants_by_sample) = match input_mode {
        InputMode::Vcf => io::vcf::read_vcf_files(&files, &criteria)?,
        InputMode::Table => io::csv::read_csv_files(&files, &criteria)?,
    };

    let rescaled = io::renormalise_saturated_sites(&mut variants_by_sample);
    if rescaled > 0 {
        eprintln!(
            "Warning: rescaled {} site(s) whose allele frequencies summed above 1 \
             (independently estimated frequencies from split multi-allelic records).",
            rescaled
        );
    }

    // FASTA fallback for reference allele (chrom-aware)
    // Only update samples that already have an entry at this position,
    // to avoid creating phantom empty entries for every sample×position.
    if let Some(ref genome) = reference {
        // A single-contig reference is unambiguous, so a contig name that does not
        // match (a table with no chrom column, or "Chromosome" vs "NC_000962.3")
        // still resolves. With several contigs there is nothing to guess from, and
        // silently reading the wrong contig's base corrupts every distance.
        let single_contig = genome.len() == 1;
        let mut unknown_chroms: std::collections::BTreeSet<&str> = std::collections::BTreeSet::new();
        let mut out_of_range = 0usize;

        for gpos in &position_set {
            let seq = match genome.get(&gpos.chrom) {
                Some(seq) => Some(seq),
                None if single_contig => genome.values().next(),
                None => {
                    unknown_chroms.insert(gpos.chrom.as_str());
                    continue;
                }
            };
            let base = seq.and_then(|s| s.get(gpos.pos.saturating_sub(1)).copied());
            let Some(b) = base else {
                out_of_range += 1;
                continue;
            };
            let ref_str = String::from(b as char);
            for sample_data in variants_by_sample.values_mut() {
                if let Some(site) = sample_data.get_mut(gpos) {
                    if site.reference_allele.is_empty() {
                        site.reference_allele = ref_str.clone();
                    }
                }
            }
        }

        if !unknown_chroms.is_empty() {
            let names: Vec<&str> = unknown_chroms.iter().copied().take(5).collect();
            eprintln!(
                "Warning: {} contig name(s) absent from the reference FASTA ({}{}). \
                 Reference alleles were not imputed there.",
                unknown_chroms.len(),
                names.join(", "),
                if unknown_chroms.len() > names.len() { ", ..." } else { "" }
            );
        }
        if out_of_range > 0 {
            eprintln!(
                "Warning: {} position(s) fall beyond the end of their contig in the reference FASTA.",
                out_of_range
            );
        }
    }

    // Check for inconsistent reference alleles across samples at the same position.
    // Collecting the distinct alleles per position, rather than comparing each one
    // against whichever was seen first, keeps both the count and the example stable:
    // the old version counted occurrences rather than positions, and which allele was
    // called "first" depended on HashMap iteration order.
    {
        use std::collections::{BTreeMap, BTreeSet};
        let mut ref_alleles_by_pos: BTreeMap<&types::GenomicPos, BTreeSet<&str>> = BTreeMap::new();
        for sample_data in variants_by_sample.values() {
            for (gpos, site) in sample_data {
                if site.reference_allele.is_empty() {
                    continue;
                }
                ref_alleles_by_pos
                    .entry(gpos)
                    .or_default()
                    .insert(site.reference_allele.as_str());
            }
        }

        let mut conflicting = ref_alleles_by_pos.iter().filter(|(_, alleles)| alleles.len() > 1);
        if let Some((gpos, alleles)) = conflicting.next() {
            let n = 1 + conflicting.count();
            let listed: Vec<&str> = alleles.iter().copied().collect();
            eprintln!(
                "Warning: {} position(s) have inconsistent REF alleles across samples, \
                 first at {}:{} ({}). \
                 Check that all inputs were called against the same reference.",
                n,
                gpos.chrom,
                gpos.pos,
                listed.join(" vs ")
            );
        }
    }

    if variants_by_sample.is_empty() {
        return Err("No samples found after filtering. Check input files and filter thresholds.".to_string());
    }

    if position_set.is_empty() {
        return Err("No polymorphic sites found after filtering.".to_string());
    }

    // Sort positions for deterministic iteration
    let all_positions: Vec<types::GenomicPos> = {
        let mut v: Vec<_> = position_set.into_iter().collect();
        v.sort();
        v
    };

    let num_loci = all_positions.len();

    eprintln!(
        "Found {} samples and {} polymorphic sites after filtering.",
        variants_by_sample.len(),
        num_loci
    );

    let samples: Vec<String> = {
        let mut s: Vec<String> = variants_by_sample.keys().cloned().collect();
        s.sort();
        s
    };
    let sample_pairs = calculation::generate_sample_pairs(&samples);

    if sample_pairs.is_empty() {
        return Err("Need at least 2 samples for pairwise distance calculation.".to_string());
    }

    let output_path_str = args.output.to_str().ok_or("Invalid output path")?;
    let formula = args.formula;
    let normalize = args.normalize;

    eprintln!(
        "Computing {:?} distances for {} pairs...",
        formula,
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

            let dist = match formula {
                Formula::Fst => calculation::fst::calculate_fst_for_pair(data1, data2, &all_positions, normalize, num_loci),
                Formula::Gst => calculation::gst::calculate_gst_for_pair(data1, data2, &all_positions),
                Formula::Nei => calculation::nei::calculate_nei_distance_for_pair(data1, data2, &all_positions),
                Formula::Chord => calculation::chord::calculate_chord_distance_for_pair(data1, data2, &all_positions, normalize, num_loci),
                Formula::BrayCurtis => calculation::bray_curtis::calculate_bray_curtis_for_pair(data1, data2, &all_positions, normalize, num_loci),
                Formula::JostD => calculation::jost_d::calculate_jost_d_for_pair(data1, data2, &all_positions, normalize, num_loci),
                Formula::Reynolds => calculation::reynolds::calculate_reynolds_distance_for_pair(data1, data2, &all_positions),
                Formula::Rogers => calculation::rogers::calculate_rogers_distance_for_pair(data1, data2, &all_positions, num_loci),
            };
            pb.inc(1);
            ((i, j), dist)
        })
        .collect();

    pb.finish_with_message("Calculation complete.");

    let distance_matrix = calculation::create_distance_matrix(&results, samples.len());
    let use_tab = output_path_str.ends_with(".tsv") || output_path_str.ends_with(".tab");

    io::csv::write_distance_matrix(output_path_str, &distance_matrix, &samples, use_tab)
        .map_err(|e| format!("Error writing output: {}", e))?;

    let elapsed = start.elapsed();
    eprintln!(
        "\n--- Summary ---\n\
         > Samples: {}\n\
         > Loci: {}\n\
         > Pairs: {}\n\
         > Formula: {:?}\n\
         > Output: {}\n\
         > Time: {:.2}s\n\
         ---------------",
        samples.len(), num_loci, sample_pairs.len(), formula, output_path_str, elapsed.as_secs_f64()
    );

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
