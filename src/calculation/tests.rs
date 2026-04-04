#[cfg(test)]
mod tests {
    use crate::calculation::common::*;
    use crate::types::{AlleleFrequencies, PositionalData, SiteData};
    use std::collections::{HashMap, HashSet};

    // ── helpers ──────────────────────────────────────────────────────────

    fn make_site(ref_allele: &str, alts: &[(&str, f64)]) -> SiteData {
        let mut freqs = AlleleFrequencies::new();
        for &(a, f) in alts {
            freqs.insert(a.to_string(), f);
        }
        SiteData {
            reference_allele: ref_allele.to_string(),
            freqs,
        }
    }

    fn make_data(entries: &[(usize, &str, &[(&str, f64)])]) -> PositionalData {
        let mut m = HashMap::new();
        for &(pos, ref_a, alts) in entries {
            m.insert(pos, make_site(ref_a, alts));
        }
        m
    }

    fn positions(poss: &[usize]) -> HashSet<usize> {
        poss.iter().copied().collect()
    }

    // ── get_all_freqs_at_pos ────────────────────────────────────────────

    #[test]
    fn freqs_ref_imputed() {
        let site = make_site("A", &[("T", 0.3)]);
        let psf = get_all_freqs_at_pos(Some(&site), None);
        assert!((psf.freqs1["A"] - 0.7).abs() < 1e-10);
        assert!((psf.freqs1["T"] - 0.3).abs() < 1e-10);
        // sample 2 has no data → ref should be 1.0
        assert!((psf.freqs2["A"] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn freqs_both_present() {
        let s1 = make_site("A", &[("T", 0.4)]);
        let s2 = make_site("A", &[("T", 0.6)]);
        let psf = get_all_freqs_at_pos(Some(&s1), Some(&s2));
        assert!((psf.freqs1["A"] - 0.6).abs() < 1e-10);
        assert!((psf.freqs2["A"] - 0.4).abs() < 1e-10);
        assert!(psf.all_alleles.contains("A"));
        assert!(psf.all_alleles.contains("T"));
    }

    #[test]
    fn freqs_clamp_negative_ref() {
        // If alt freqs sum > 1 (rounding), ref should be 0 not negative
        let site = make_site("A", &[("T", 0.7), ("G", 0.4)]);
        let psf = get_all_freqs_at_pos(Some(&site), None);
        assert!(psf.freqs1["A"] >= 0.0);
    }

    // ── heterozygosity / homozygosity ───────────────────────────────────

    #[test]
    fn het_monomorphic() {
        let mut f = AlleleFrequencies::new();
        f.insert("A".into(), 1.0);
        assert!((heterozygosity(&f) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn het_biallelic_equal() {
        let mut f = AlleleFrequencies::new();
        f.insert("A".into(), 0.5);
        f.insert("T".into(), 0.5);
        assert!((heterozygosity(&f) - 0.5).abs() < 1e-10);
    }

    // ── FST ─────────────────────────────────────────────────────────────

    #[test]
    fn fst_identical_is_zero() {
        let d1 = make_data(&[(100, "A", &[("T", 0.3)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.3)])]);
        let pos = positions(&[100]);
        let fst = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &pos, false, 1);
        assert!(fst.abs() < 1e-10, "Identical pops should have FST=0, got {}", fst);
    }

    #[test]
    fn fst_fixed_difference_is_one() {
        // Pop1: 100% T, Pop2: 100% A (ref)
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        let pos = positions(&[100]);
        let fst = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &pos, false, 1);
        assert!((fst - 1.0).abs() < 1e-10, "Fixed difference should give FST=1, got {}", fst);
    }

    #[test]
    fn fst_normalized() {
        let d1 = make_data(&[(100, "A", &[("T", 0.8)]), (200, "C", &[("G", 0.8)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.2)]), (200, "C", &[("G", 0.2)])]);
        let pos = positions(&[100, 200]);
        let raw = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &pos, false, 2);
        let norm = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &pos, true, 2);
        assert!((norm - raw / 2.0).abs() < 1e-10);
    }

    // ── GST ─────────────────────────────────────────────────────────────

    #[test]
    fn gst_identical_is_zero() {
        let d1 = make_data(&[(100, "A", &[("T", 0.5)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.5)])]);
        let pos = positions(&[100]);
        let gst = crate::calculation::gst::calculate_gst_for_pair(&d1, &d2, &pos);
        assert!(gst.abs() < 1e-10);
    }

    // ── Jost's D ────────────────────────────────────────────────────────

    #[test]
    fn jost_d_identical_is_zero() {
        let d1 = make_data(&[(100, "A", &[("T", 0.3)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.3)])]);
        let pos = positions(&[100]);
        let d = crate::calculation::jost_d::calculate_jost_d_for_pair(&d1, &d2, &pos, false, 1);
        assert!(d.abs() < 1e-10, "Identical → D=0, got {}", d);
    }

    #[test]
    fn jost_d_fixed_difference() {
        // Pop1: 100% alt, Pop2: 100% ref → max differentiation
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        let pos = positions(&[100]);
        let d = crate::calculation::jost_d::calculate_jost_d_for_pair(&d1, &d2, &pos, false, 1);
        // For 2 pops with fixed alleles: Hs=0, Ht=0.5 → D = 2*(0.5-0)/(1-0) = 1.0
        assert!((d - 1.0).abs() < 1e-10, "Fixed diff → D=1, got {}", d);
    }

    #[test]
    fn jost_d_is_always_non_negative() {
        // Any freq combo should give D >= 0
        for af1 in &[0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0] {
            for af2 in &[0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0] {
                let d1 = make_data(&[(100, "A", &[("T", *af1)])]);
                let d2 = make_data(&[(100, "A", &[("T", *af2)])]);
                let pos = positions(&[100]);
                let d = crate::calculation::jost_d::calculate_jost_d_for_pair(&d1, &d2, &pos, false, 1);
                assert!(d >= -1e-10, "Jost D should be >= 0, got {} for af1={} af2={}", d, af1, af2);
            }
        }
    }

    // ── Nei's D ─────────────────────────────────────────────────────────

    #[test]
    fn nei_identical_is_zero() {
        let d1 = make_data(&[(100, "A", &[("T", 0.4)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.4)])]);
        let pos = positions(&[100]);
        let d = crate::calculation::nei::calculate_nei_distance_for_pair(&d1, &d2, &pos);
        assert!(d.abs() < 1e-10, "Identical → D=0, got {}", d);
    }

    #[test]
    fn nei_fixed_diff_is_inf() {
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        let pos = positions(&[100]);
        let d = crate::calculation::nei::calculate_nei_distance_for_pair(&d1, &d2, &pos);
        assert!(d.is_infinite(), "Fixed diff → D=inf, got {}", d);
    }

    // ── Chord ───────────────────────────────────────────────────────────

    #[test]
    fn chord_identical_is_zero() {
        let d1 = make_data(&[(100, "A", &[("T", 0.5)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.5)])]);
        let pos = positions(&[100]);
        let d = crate::calculation::chord::calculate_chord_distance_for_pair(&d1, &d2, &pos, false, 1);
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn chord_fixed_diff() {
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        let pos = positions(&[100]);
        let d = crate::calculation::chord::calculate_chord_distance_for_pair(&d1, &d2, &pos, false, 1);
        // sqrt(2 * (1 - 0)) = sqrt(2) ≈ 1.414
        assert!((d - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    // ── Bray-Curtis ─────────────────────────────────────────────────────

    #[test]
    fn bray_curtis_identical_is_zero() {
        let d1 = make_data(&[(100, "A", &[("T", 0.3)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.3)])]);
        let pos = positions(&[100]);
        let d = crate::calculation::bray_curtis::calculate_bray_curtis_for_pair(&d1, &d2, &pos, false, 1);
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn bray_curtis_fixed_diff() {
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        let pos = positions(&[100]);
        let d = crate::calculation::bray_curtis::calculate_bray_curtis_for_pair(&d1, &d2, &pos, false, 1);
        // 0.5 * (|1-0| + |0-1|) = 0.5 * 2 = 1.0
        assert!((d - 1.0).abs() < 1e-10);
    }

    // ── Reynolds ────────────────────────────────────────────────────────

    #[test]
    fn reynolds_identical_is_zero() {
        let d1 = make_data(&[(100, "A", &[("T", 0.5)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.5)])]);
        let pos = positions(&[100]);
        let d = crate::calculation::reynolds::calculate_reynolds_distance_for_pair(&d1, &d2, &pos);
        assert!(d.abs() < 1e-10);
    }

    // ── Rogers ──────────────────────────────────────────────────────────

    #[test]
    fn rogers_identical_is_zero() {
        let d1 = make_data(&[(100, "A", &[("T", 0.5)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.5)])]);
        let pos = positions(&[100]);
        let d = crate::calculation::rogers::calculate_rogers_distance_for_pair(&d1, &d2, &pos, 1);
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn rogers_fixed_diff() {
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        let pos = positions(&[100]);
        let d = crate::calculation::rogers::calculate_rogers_distance_for_pair(&d1, &d2, &pos, 1);
        // sqrt((1+1) / 2) = 1.0
        assert!((d - 1.0).abs() < 1e-10);
    }

    // ── generate_sample_pairs ───────────────────────────────────────────

    #[test]
    fn pairs_count() {
        let samples: Vec<String> = (0..5).map(|i| format!("s{}", i)).collect();
        let pairs = crate::calculation::generate_sample_pairs(&samples);
        assert_eq!(pairs.len(), 10); // 5 choose 2
    }

    #[test]
    fn pairs_empty() {
        let samples: Vec<String> = vec![];
        let pairs = crate::calculation::generate_sample_pairs(&samples);
        assert!(pairs.is_empty());
    }

    // ── distance matrix ─────────────────────────────────────────────────

    #[test]
    fn matrix_symmetric() {
        let results = vec![((0, 1), 0.5), ((0, 2), 0.3), ((1, 2), 0.7)];
        let m = crate::calculation::create_distance_matrix(&results, 3);
        assert_eq!(m[0][1], m[1][0]);
        assert_eq!(m[0][2], m[2][0]);
        assert_eq!(m[1][2], m[2][1]);
        assert_eq!(m[0][0], 0.0);
    }

    // ── Multi-locus consistency ─────────────────────────────────────────

    #[test]
    fn fst_multi_locus_sum() {
        // 2 loci with known per-site FST
        let d1 = make_data(&[
            (100, "A", &[("T", 0.9)]),
            (200, "C", &[("G", 0.1)]),
        ]);
        let d2 = make_data(&[
            (100, "A", &[("T", 0.1)]),
            (200, "C", &[("G", 0.9)]),
        ]);
        let pos = positions(&[100, 200]);
        let fst = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &pos, false, 2);
        // Each site should contribute equally → sum should be 2x single site
        let single_pos = positions(&[100]);
        let fst_single = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &single_pos, false, 1);
        assert!((fst - 2.0 * fst_single).abs() < 1e-10);
    }
}
