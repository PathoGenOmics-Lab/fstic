use crate::types::{AlleleFrequencies, PositionalData, SiteData};
use rayon::prelude::*;
use std::collections::{HashSet};

/// Calculates Nei's standard genetic distance D = -ln(I).
/// This is a global metric and is not affected by the --normalize flag.
pub fn calculate_nei_distance_for_pair(
    data1: &PositionalData,
    data2: &PositionalData,
    all_positions: &HashSet<usize>,
) -> f64 {
    // Sum the components of genetic identity over all loci in parallel
    let (j_xy, j_x, j_y) = all_positions
        .par_iter()
        .map(|&pos| {
            let (all_alleles, freqs1, freqs2) = get_all_freqs_at_pos(data1.get(&pos), data2.get(&pos));

            let mut site_j_xy = 0.0;
            let mut site_j_x = 0.0;
            let mut site_j_y = 0.0;

            for allele in all_alleles {
                let p_i = *freqs1.get(&allele).unwrap_or(&0.0);
                let q_i = *freqs2.get(&allele).unwrap_or(&0.0);
                site_j_xy += p_i * q_i;
                site_j_x += p_i.powi(2);
                site_j_y += q_i.powi(2);
            }
            (site_j_xy, site_j_x, site_j_y)
        })
        .reduce(|| (0.0, 0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2));

    // Calculate Nei's genetic identity I
    let denominator = (j_x * j_y).sqrt();
    if denominator > 0.0 {
        let identity = j_xy / denominator;
        if identity > 0.0 {
            -identity.ln()
        } else {
            f64::INFINITY // Distance is infinite if identity is zero
        }
    } else {
        f64::INFINITY // Or if one population has no variation
    }
}

/// Helper to get complete frequency maps for both samples at a site.
fn get_all_freqs_at_pos<'a>(
    site1_data: Option<&'a SiteData>,
    site2_data: Option<&'a SiteData>,
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
