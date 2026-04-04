use crate::types::{AlleleFrequencies, SiteData};
use std::collections::HashSet;

/// Get complete frequency maps for both samples at a site.
/// Returns the set of all alleles and the frequency maps for each sample,
/// with the reference allele frequency computed as 1 - sum(alt freqs).
pub fn get_all_freqs_at_pos(
    site1_data: Option<&SiteData>,
    site2_data: Option<&SiteData>,
) -> (HashSet<String>, AlleleFrequencies, AlleleFrequencies) {
    let mut all_alleles: HashSet<String> = HashSet::new();
    let mut freqs1 = AlleleFrequencies::new();
    let mut freqs2 = AlleleFrequencies::new();

    let mut ref_allele = String::new();
    if let Some(data) = site1_data {
        ref_allele = data.reference_allele.clone();
    } else if let Some(data) = site2_data {
        ref_allele = data.reference_allele.clone();
    }

    if let Some(data) = site1_data {
        freqs1 = data.freqs.clone();
        all_alleles.extend(data.freqs.keys().cloned());
    }

    if let Some(data) = site2_data {
        freqs2 = data.freqs.clone();
        all_alleles.extend(data.freqs.keys().cloned());
    }

    if !ref_allele.is_empty() {
        all_alleles.insert(ref_allele.clone());
    }

    let sum1: f64 = freqs1.values().sum();
    if !ref_allele.is_empty() {
        freqs1.insert(ref_allele.clone(), 1.0 - sum1);
    }

    let sum2: f64 = freqs2.values().sum();
    if !ref_allele.is_empty() {
        freqs2.insert(ref_allele, 1.0 - sum2);
    }

    (all_alleles, freqs1, freqs2)
}
