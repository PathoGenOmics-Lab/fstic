use crate::types::{FilterCriteria, SampleVariants, VcfVariant};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

pub fn read_vcf_files(
    files: &[PathBuf],
    criteria: &FilterCriteria,
) -> (HashSet<usize>, SampleVariants) {
    let multiallelic_count = AtomicUsize::new(0);

    let all_variants: Vec<VcfVariant> = files
        .par_iter()
        .flat_map(|file| parse_and_filter_vcf(file, criteria, &multiallelic_count))
        .collect();

    let skipped = multiallelic_count.load(Ordering::Relaxed);
    if skipped > 0 {
        eprintln!(
            "Warning: {} multi-allelic site(s) had ALT alleles without matching per-allele FREQ values; \
             those alleles were assigned frequency 0.",
            skipped
        );
    }

    let mut variants_by_sample: SampleVariants = HashMap::new();
    let mut all_positions: HashSet<usize> = HashSet::new();

    for variant in all_variants {
        all_positions.insert(variant.pos);
        let sample_map = variants_by_sample
            .entry(variant.sample.clone())
            .or_default();

        let site_data = sample_map
            .entry(variant.pos)
            .or_default();

        site_data.reference_allele = variant.ref_allele;
        site_data.freqs.insert(variant.alt_allele, variant.alt_freq.unwrap_or(0.0));
    }

    (all_positions, variants_by_sample)
}

fn parse_and_filter_vcf(file: &Path, criteria: &FilterCriteria, multiallelic_count: &AtomicUsize) -> Vec<VcfVariant> {
    let sample_id = file
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown_sample")
        .to_string();

    let vcf_file = File::open(file).expect("Cannot open VCF file");
    let reader = BufReader::new(vcf_file);
    let warnings = Arc::new(Mutex::new(HashSet::new()));

    let mut variants = Vec::new();
    for line in reader.lines().flatten() {
        if line.starts_with('#') { continue; }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 { continue; }

        let parsed = parse_vcf_line(&fields, &sample_id, file, &warnings, multiallelic_count);
        variants.extend(parsed);
    }

    variants.into_iter().filter(|v| {
        let depth_ok = v.total_dp.map_or(true, |dp| dp >= criteria.min_depth);
        let freq_ok = v.alt_freq.map_or(true, |freq| freq >= criteria.min_freq);
        let alt_reads_ok = v.alt_dp.map_or(true, |ad| ad >= criteria.min_alt_reads);
        let alt_rev_reads_ok = v.alt_rv.map_or(true, |arv| arv >= criteria.min_alt_rev_reads);
        depth_ok && freq_ok && alt_reads_ok && alt_rev_reads_ok
    }).collect()
}

fn parse_vcf_line(
    fields: &[&str],
    sample_id: &str,
    filepath: &Path,
    warnings: &Arc<Mutex<HashSet<String>>>,
    multiallelic_count: &AtomicUsize,
) -> Vec<VcfVariant> {
    let pos = match fields[1].parse::<usize>() {
        Ok(p) => p,
        Err(_) => return vec![],
    };
    let ref_allele = fields[3].to_string();
    let alt_alleles: Vec<&str> = fields[4].split(',').collect();

    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let sample_values: Vec<&str> = fields[9].split(':').collect();
    let format_map: HashMap<&str, &str> = format_keys.iter().copied().zip(sample_values.iter().copied()).collect();

    let get_value = |key: &str| -> &str {
        format_map.get(key).copied().unwrap_or(".")
    };

    let parse_and_warn = |key: &str| -> Option<u32> {
        match get_value(key) {
            "." => {
                let mut locked_warnings = warnings.lock().unwrap();
                if !locked_warnings.contains(key) {
                    eprintln!("Warning: VCF file '{}' is missing FORMAT field '{}'. Filtering on this field will be skipped.", filepath.display(), key);
                    locked_warnings.insert(key.to_string());
                }
                None
            },
            val => val.parse::<u32>().ok(),
        }
    };

    let total_dp = parse_and_warn("DP");

    // Parse per-allele AD values (comma-separated: ref,alt1,alt2,...)
    let ad_values: Vec<Option<u32>> = match get_value("AD") {
        "." => {
            let mut locked_warnings = warnings.lock().unwrap();
            if !locked_warnings.contains("AD") {
                eprintln!("Warning: VCF file '{}' is missing FORMAT field 'AD'. Filtering on this field will be skipped.", filepath.display());
                locked_warnings.insert("AD".to_string());
            }
            vec![None; alt_alleles.len()]
        },
        val => {
            let parts: Vec<Option<u32>> = val.split(',').map(|v| v.parse::<u32>().ok()).collect();
            // Skip first element (ref AD), take alt ADs
            if parts.len() > 1 {
                parts[1..].to_vec()
            } else {
                // Single value: apply to first alt allele
                vec![parts.first().copied().flatten(); alt_alleles.len()]
            }
        }
    };

    let alt_rv = parse_and_warn("ADR");

    // Parse per-allele FREQ values (comma-separated)
    let freq_raw = get_value("FREQ");
    let freq_values: Vec<Option<f64>> = match freq_raw {
        "." => {
            let mut locked_warnings = warnings.lock().unwrap();
            if !locked_warnings.contains("FREQ") {
                eprintln!("Warning: VCF file '{}' is missing FORMAT field 'FREQ'.", filepath.display());
                locked_warnings.insert("FREQ".to_string());
            }
            vec![None; alt_alleles.len()]
        },
        val => {
            val.split(',').map(|v| {
                let has_percent = v.contains('%');
                v.trim_end_matches('%').parse::<f64>().ok().map(|f| {
                    if has_percent { f / 100.0 } else { f }
                })
            }).collect()
        }
    };

    let is_multiallelic = alt_alleles.len() > 1;
    let mut variants = Vec::with_capacity(alt_alleles.len());

    for (i, alt) in alt_alleles.iter().enumerate() {
        let alt_dp_i = ad_values.get(i).copied().flatten();
        let alt_freq_i = freq_values.get(i).copied().flatten();

        if is_multiallelic && alt_freq_i.is_none() {
            multiallelic_count.fetch_add(1, Ordering::Relaxed);
        }

        variants.push(VcfVariant {
            sample: sample_id.to_string(),
            pos,
            ref_allele: ref_allele.clone(),
            alt_allele: alt.to_string(),
            total_dp,
            alt_dp: alt_dp_i,
            alt_freq: alt_freq_i,
            alt_rv,
        });
    }

    variants
}
