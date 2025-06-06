use crate::types::{FilterCriteria, SampleVariants, SiteData, VcfVariant};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

pub fn read_vcf_files(
    files: &[PathBuf],
    criteria: &FilterCriteria,
) -> (HashSet<usize>, SampleVariants) {
    let all_variants: Vec<VcfVariant> = files
        .par_iter()
        .flat_map(|file| parse_and_filter_vcf(file, criteria))
        .collect();

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

fn parse_and_filter_vcf(file: &Path, criteria: &FilterCriteria) -> Vec<VcfVariant> {
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

        if let Some(variant) = parse_vcf_line(&fields, &sample_id, file, &warnings) {
            variants.push(variant);
        }
    }

    variants.into_iter().filter(|v| {
        let depth_ok = v.total_dp.map_or(true, |dp| dp >= criteria.min_depth);
        let freq_ok = v.alt_freq.map_or(true, |freq| freq >= criteria.min_freq);
        let alt_reads_ok = v.alt_dp.map_or(true, |ad| ad >= criteria.min_alt_reads);
        let alt_rev_reads_ok = v.alt_rv.map_or(true, |arv| arv >= criteria.min_alt_rev_reads);
        depth_ok && freq_ok && alt_reads_ok && alt_rev_reads_ok
    }).collect()
}

fn parse_vcf_line<'a>(
    fields: &[&'a str],
    sample_id: &str,
    filepath: &Path,
    warnings: &Arc<Mutex<HashSet<String>>>,
) -> Option<VcfVariant> {
    let pos = fields[1].parse::<usize>().ok()?;
    let ref_allele = fields[3].to_string();
    let alt_alleles = fields[4].split(',');

    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let sample_values: Vec<&str> = fields[9].split(':').collect();
    let format_map: HashMap<&str, &str> = format_keys.iter().copied().zip(sample_values.iter().copied()).collect();

    let get_value = |key: &str| *format_map.get(key).unwrap_or(&".");

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
    let alt_dp = parse_and_warn("AD");
    let alt_rv = parse_and_warn("ADR");

    let alt_freq = match get_value("FREQ") {
        "." => {
            let mut locked_warnings = warnings.lock().unwrap();
            if !locked_warnings.contains("FREQ") {
                 eprintln!("Warning: VCF file '{}' is missing FORMAT field 'FREQ'.", filepath.display());
                 locked_warnings.insert("FREQ".to_string());
            }
            None
        },
        val => val.trim_end_matches('%').parse::<f64>().ok().map(|f| f / 100.0)
    };

    if alt_alleles.clone().count() == 1 {
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
    } else {
        None
    }
}
