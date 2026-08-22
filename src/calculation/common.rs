use crate::types::{AlleleFrequencies, GenomicPos, SiteData};
use rayon::prelude::*;
use std::collections::HashSet;

/// Loci per reduction chunk.
///
/// rayon splits an indexed iterator according to the current thread count and work
/// stealing, so a plain `par_iter().sum()` adds the per-locus terms in an order that
/// varies with `--workers`. Floating point addition is not associative, so the same
/// input gave different output bytes on different machines: on 320k loci, `-w 3`
/// produced 91356.2127736697 where every other thread count gave ...699.
///
/// A fixed chunk size makes the summation tree a function of the input length alone.
const REDUCTION_CHUNK: usize = 4096;

/// Sums a per-locus quantity in an order that does not depend on the thread count.
pub fn sum_per_locus<F>(all_positions: &[GenomicPos], per_locus: F) -> f64
where
    F: Fn(&GenomicPos) -> f64 + Sync + Send,
{
    all_positions
        .par_chunks(REDUCTION_CHUNK)
        .map(|chunk| chunk.iter().map(&per_locus).sum::<f64>())
        .collect::<Vec<f64>>()
        .into_iter()
        .sum()
}

/// Same as [`sum_per_locus`] for metrics that accumulate several running totals.
pub fn sum_per_locus_n<const N: usize, F>(all_positions: &[GenomicPos], per_locus: F) -> [f64; N]
where
    F: Fn(&GenomicPos) -> [f64; N] + Sync + Send,
{
    all_positions
        .par_chunks(REDUCTION_CHUNK)
        .map(|chunk| {
            let mut acc = [0.0_f64; N];
            for pos in chunk {
                let terms = per_locus(pos);
                for (a, t) in acc.iter_mut().zip(terms) {
                    *a += t;
                }
            }
            acc
        })
        .collect::<Vec<[f64; N]>>()
        .into_iter()
        .fold([0.0_f64; N], |mut acc, chunk| {
            for (a, t) in acc.iter_mut().zip(chunk) {
                *a += t;
            }
            acc
        })
}

/// Stand-in for the reference base at a site whose REF is unknown.
///
/// The position list is the union of variant sites over *all* samples, so for any
/// given pair there are usually positions where neither member has a record. Both
/// samples are reference there, but no REF allele was ever read for them. Angle
/// brackets cannot occur in a called SNP allele, so this label never collides with
/// real data.
const UNKNOWN_REF: &str = "<REF>";

/// Per-site frequency data for a pair of samples, with reference allele imputed.
pub struct PairSiteFreqs {
    pub all_alleles: HashSet<String>,
    pub freqs1: AlleleFrequencies,
    pub freqs2: AlleleFrequencies,
}

/// Builds complete allele frequency profiles for two samples at a single site.
///
/// For each sample the reference allele frequency is imputed as `1 - sum(alt_freqs)`,
/// so a sample with no record at the site is treated as homozygous reference.
/// Returns the union of all alleles observed plus both frequency maps.
pub fn get_all_freqs_at_pos(
    site1_data: Option<&SiteData>,
    site2_data: Option<&SiteData>,
) -> PairSiteFreqs {
    let mut all_alleles: HashSet<String> = HashSet::new();
    let mut freqs1 = AlleleFrequencies::new();
    let mut freqs2 = AlleleFrequencies::new();

    // Determine the reference allele from whichever sample has it. Neither sample
    // having a record here means the site is private to some other sample: both are
    // reference, so a placeholder keeps the profiles complete and the site neutral.
    let ref_allele = site1_data
        .map(|d| &d.reference_allele)
        .filter(|r| !r.is_empty())
        .or_else(|| {
            site2_data
                .map(|d| &d.reference_allele)
                .filter(|r| !r.is_empty())
        })
        .cloned()
        .unwrap_or_else(|| UNKNOWN_REF.to_string());

    if let Some(data) = site1_data {
        freqs1.clone_from(&data.freqs);
        all_alleles.extend(data.freqs.keys().cloned());
    }

    if let Some(data) = site2_data {
        freqs2.clone_from(&data.freqs);
        all_alleles.extend(data.freqs.keys().cloned());
    }

    all_alleles.insert(ref_allele.clone());

    // Only impute ref freq if the sample doesn't already have an explicit
    // frequency for the ref allele (e.g. from table input where ref is
    // listed as an allele with its own frequency).
    if !freqs1.contains_key(&ref_allele) {
        let sum1: f64 = freqs1.values().sum();
        freqs1.insert(ref_allele.clone(), (1.0 - sum1).max(0.0));
    }

    if !freqs2.contains_key(&ref_allele) {
        let sum2: f64 = freqs2.values().sum();
        freqs2.insert(ref_allele, (1.0 - sum2).max(0.0));
    }

    PairSiteFreqs {
        all_alleles,
        freqs1,
        freqs2,
    }
}

/// Heterozygosity (gene diversity) for a frequency map: `1 - sum(p_i^2)`.
#[inline]
pub fn heterozygosity(freqs: &AlleleFrequencies) -> f64 {
    1.0 - freqs.values().map(|p| p.powi(2)).sum::<f64>()
}

/// Homozygosity for a frequency map: `sum(p_i^2)`.
#[inline]
#[allow(dead_code)]
pub fn homozygosity(freqs: &AlleleFrequencies) -> f64 {
    freqs.values().map(|p| p.powi(2)).sum::<f64>()
}

/// Pooled-population heterozygosity for two samples at one site.
pub fn pooled_heterozygosity(psf: &PairSiteFreqs) -> f64 {
    1.0 - psf
        .all_alleles
        .iter()
        .map(|allele| {
            let p = psf.freqs1.get(allele).copied().unwrap_or(0.0);
            let q = psf.freqs2.get(allele).copied().unwrap_or(0.0);
            ((p + q) / 2.0).powi(2)
        })
        .sum::<f64>()
}
