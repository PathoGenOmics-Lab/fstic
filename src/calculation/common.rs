use crate::types::{AlleleFrequencies, SiteData};
use std::collections::HashSet;

/// Per-site frequency data for a pair of samples, with reference allele imputed.
pub struct PairSiteFreqs {
    pub all_alleles: HashSet<String>,
    pub freqs1: AlleleFrequencies,
    pub freqs2: AlleleFrequencies,
}

/// Builds complete allele frequency profiles for two samples at a single site.
///
/// For each sample the reference allele frequency is imputed as `1 - sum(alt_freqs)`.
/// Returns the union of all alleles observed plus both frequency maps.
pub fn get_all_freqs_at_pos(
    site1_data: Option<&SiteData>,
    site2_data: Option<&SiteData>,
) -> PairSiteFreqs {
    let mut all_alleles: HashSet<String> = HashSet::new();
    let mut freqs1 = AlleleFrequencies::new();
    let mut freqs2 = AlleleFrequencies::new();

    // Determine the reference allele from whichever sample has it
    let ref_allele = site1_data
        .map(|d| &d.reference_allele)
        .filter(|r| !r.is_empty())
        .or_else(|| site2_data.map(|d| &d.reference_allele).filter(|r| !r.is_empty()))
        .cloned()
        .unwrap_or_default();

    if let Some(data) = site1_data {
        freqs1.clone_from(&data.freqs);
        all_alleles.extend(data.freqs.keys().cloned());
    }

    if let Some(data) = site2_data {
        freqs2.clone_from(&data.freqs);
        all_alleles.extend(data.freqs.keys().cloned());
    }

    if !ref_allele.is_empty() {
        all_alleles.insert(ref_allele.clone());

        let sum1: f64 = freqs1.values().sum();
        freqs1.insert(ref_allele.clone(), (1.0 - sum1).max(0.0));

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
