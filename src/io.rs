use crate::types::SampleVariants;
use std::fs::File;
use std::io::{BufRead, BufReader, Error};
use std::path::PathBuf;

pub mod csv;
pub mod fasta;
pub mod vcf;

/// Rescales any site whose allele frequencies sum above 1 back to a distribution,
/// returning how many sites were touched.
///
/// Split multi-allelic records (what `bcftools norm -m-any` produces) carry
/// independently estimated frequencies that can add up to more than 1. The reference
/// imputation then clamps to 0, which hides the problem without solving it: the
/// vector is no longer a distribution, `1 - sum(p^2)` can go negative, and every
/// metric leaves its theoretical range. One locus like that was enough to push FST to
/// 1.64 and Rogers to 1.06.
///
/// The tolerance keeps ordinary rounding noise silent.
pub fn renormalise_saturated_sites(variants_by_sample: &mut SampleVariants) -> usize {
    const TOLERANCE: f64 = 1e-6;
    let mut rescaled = 0;

    for sample_data in variants_by_sample.values_mut() {
        for site in sample_data.values_mut() {
            let sum: f64 = site.freqs.values().sum();
            if sum > 1.0 + TOLERANCE {
                for freq in site.freqs.values_mut() {
                    *freq /= sum;
                }
                rescaled += 1;
            }
        }
    }

    rescaled
}

pub fn read_file_list(path: &str) -> Result<Vec<PathBuf>, Error> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut paths = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if !trimmed.is_empty() && !trimmed.starts_with('#') {
            paths.push(PathBuf::from(trimmed));
        }
    }
    Ok(paths)
}
