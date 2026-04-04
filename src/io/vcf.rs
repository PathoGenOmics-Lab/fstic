use crate::types::{FilterCriteria, SampleVariants, VcfVariant};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

pub fn read_vcf_files(
    files: &[PathBuf],
    criteria: &FilterCriteria,
) -> (HashSet<usize>, SampleVariants) {
    let multi_allelic_count = AtomicUsize::new(0);

    let all_variants: Vec<VcfVariant> = files
        .par_iter()
        .flat_map(|file| parse_and_filter_vcf(file, criteria, &multi_allelic_count))
        .collect();

    let skipped = multi_allelic_count.load(Ordering::Relaxed);
    if skipped > 0 {
        eprintln!(
            "Warning: skipped {} multi-allelic site(s) (>1 ALT allele). \
             Only bi-allelic sites are currently supported.",
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

        let site_data = sample_map.entry(variant.pos).or_default();
        site_data.reference_allele = variant.ref_allele;
        site_data
            .freqs
            .insert(variant.alt_allele, variant.alt_freq.unwrap_or(0.0));
    }

    (all_positions, variants_by_sample)
}

fn parse_and_filter_vcf(
    file: &Path,
    criteria: &FilterCriteria,
    multi_allelic_count: &AtomicUsize,
) -> Vec<VcfVariant> {
    let sample_id = file
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown_sample")
        .to_string();

    let vcf_file = match File::open(file) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error: cannot open VCF file {}: {}", file.display(), e);
            return Vec::new();
        }
    };
    let reader = BufReader::new(vcf_file);
    let mut warned: HashSet<String> = HashSet::new();

    let mut variants = Vec::new();
    for line in reader.lines().map_while(Result::ok) {
        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            continue;
        }

        // Skip multi-allelic sites with a count
        if fields[4].contains(',') {
            multi_allelic_count.fetch_add(1, Ordering::Relaxed);
            continue;
        }

        if let Some(variant) = parse_vcf_line(&fields, &sample_id, file, &mut warned) {
            variants.push(variant);
        }
    }

    variants
        .into_iter()
        .filter(|v| {
            let depth_ok = v.total_dp.is_none_or(|dp| dp >= criteria.min_depth);
            let freq_ok = v.alt_freq.is_none_or(|freq| freq >= criteria.min_freq);
            let alt_reads_ok = v.alt_dp.is_none_or(|ad| ad >= criteria.min_alt_reads);
            let alt_rev_reads_ok = v
                .alt_rv
                .is_none_or(|arv| arv >= criteria.min_alt_rev_reads);
            depth_ok && freq_ok && alt_reads_ok && alt_rev_reads_ok
        })
        .collect()
}

fn parse_vcf_line(
    fields: &[&str],
    sample_id: &str,
    filepath: &Path,
    warned: &mut HashSet<String>,
) -> Option<VcfVariant> {
    let pos = fields[1].parse::<usize>().ok()?;
    let ref_allele = fields[3].to_string();

    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let sample_values: Vec<&str> = fields[9].split(':').collect();
    let format_map: HashMap<&str, &str> = format_keys
        .iter()
        .copied()
        .zip(sample_values.iter().copied())
        .collect();

    let get_value = |key: &str| -> &str { format_map.get(key).copied().unwrap_or(".") };

    let parse_and_warn = |key: &str, warned: &mut HashSet<String>| -> Option<u32> {
        match get_value(key) {
            "." => {
                if warned.insert(key.to_string()) {
                    eprintln!(
                        "Warning: VCF file '{}' is missing FORMAT field '{}'. \
                         Filtering on this field will be skipped.",
                        filepath.display(),
                        key
                    );
                }
                None
            }
            val => val.parse::<u32>().ok(),
        }
    };

    let total_dp = parse_and_warn("DP", warned);
    let alt_dp = parse_and_warn("AD", warned);
    let alt_rv = parse_and_warn("ADR", warned);

    // FIX: only divide by 100 when the raw value actually contains '%'
    let alt_freq = match get_value("FREQ") {
        "." => {
            if warned.insert("FREQ".to_string()) {
                eprintln!(
                    "Warning: VCF file '{}' is missing FORMAT field 'FREQ'.",
                    filepath.display()
                );
            }
            None
        }
        val => {
            if val.ends_with('%') {
                val.trim_end_matches('%')
                    .parse::<f64>()
                    .ok()
                    .map(|f| f / 100.0)
            } else {
                val.parse::<f64>().ok()
            }
        }
    };

    Some(VcfVariant {
        sample: sample_id.to_string(),
        pos,
        ref_allele,
        alt_allele: fields[4].to_string(),
        total_dp,
        alt_dp,
        alt_freq,
        alt_rv,
    })
}
