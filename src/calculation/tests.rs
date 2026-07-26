#[cfg(test)]
mod tests {
    use crate::calculation::common::*;
    use crate::types::{AlleleFrequencies, GenomicPos, PositionalData, SiteData};
    use std::collections::HashMap;

    fn gp(pos: usize) -> GenomicPos {
        GenomicPos { chrom: "chr1".into(), pos }
    }

    fn make_site(ref_allele: &str, alts: &[(&str, f64)]) -> SiteData {
        let mut freqs = AlleleFrequencies::new();
        for &(a, f) in alts {
            freqs.insert(a.to_string(), f);
        }
        SiteData { reference_allele: ref_allele.to_string(), freqs }
    }

    fn make_data(entries: &[(usize, &str, &[(&str, f64)])]) -> PositionalData {
        let mut m = HashMap::new();
        for &(pos, ref_a, alts) in entries {
            m.insert(gp(pos), make_site(ref_a, alts));
        }
        m
    }

    fn gps(poss: &[usize]) -> Vec<GenomicPos> {
        poss.iter().map(|&p| gp(p)).collect()
    }

    // ── get_all_freqs_at_pos ────────────────────────────────────────────

    #[test]
    fn freqs_ref_imputed() {
        let site = make_site("A", &[("T", 0.3)]);
        let psf = get_all_freqs_at_pos(Some(&site), None);
        assert!((psf.freqs1["A"] - 0.7).abs() < 1e-10);
        assert!((psf.freqs1["T"] - 0.3).abs() < 1e-10);
        assert!((psf.freqs2["A"] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn freqs_both_present() {
        let s1 = make_site("A", &[("T", 0.4)]);
        let s2 = make_site("A", &[("T", 0.6)]);
        let psf = get_all_freqs_at_pos(Some(&s1), Some(&s2));
        assert!((psf.freqs1["A"] - 0.6).abs() < 1e-10);
        assert!((psf.freqs2["A"] - 0.4).abs() < 1e-10);
    }

    #[test]
    fn freqs_clamp_negative_ref() {
        let site = make_site("A", &[("T", 0.7), ("G", 0.4)]);
        let psf = get_all_freqs_at_pos(Some(&site), None);
        assert!(psf.freqs1["A"] >= 0.0);
    }

    #[test]
    fn freqs_no_ref() {
        // REF unknown: the missing mass still has to go somewhere, so a placeholder
        // reference allele carries it and the profile stays a distribution.
        let s1 = SiteData { reference_allele: String::new(), freqs: { let mut f = AlleleFrequencies::new(); f.insert("T".into(), 0.5); f } };
        let psf = get_all_freqs_at_pos(Some(&s1), None);
        assert_eq!(psf.freqs1.len(), 2);
        let sum: f64 = psf.freqs1.values().sum();
        assert!((sum - 1.0).abs() < 1e-10);
    }

    #[test]
    fn freqs_both_none() {
        // Site private to some other sample: both members of this pair are reference.
        let psf = get_all_freqs_at_pos(None, None);
        assert_eq!(psf.all_alleles.len(), 1);
        let a = psf.all_alleles.iter().next().unwrap();
        assert!((psf.freqs1[a] - 1.0).abs() < 1e-10);
        assert!((psf.freqs2[a] - 1.0).abs() < 1e-10);
    }

    /// A locus where neither sample of the pair has a record is a locus where both are
    /// reference, so it must add nothing to the cumulative metrics. Without this, a
    /// pairwise distance changes when an unrelated sample joins the run.
    ///
    /// Nei's D and Rogers are excluded on purpose: both are defined over the whole
    /// locus set (a ratio and a per-locus mean respectively), so a shared monomorphic
    /// locus legitimately moves them. They are covered by
    /// `identical_samples_zero_despite_foreign_loci` instead.
    #[test]
    fn absent_site_is_neutral_for_cumulative_metrics() {
        let d1 = make_data(&[(100, "A", &[("T", 0.7)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.2)])]);
        let one = gps(&[100]);
        // Position 999 is private to a third sample and absent from both d1 and d2.
        let with_absent = gps(&[100, 999]);

        macro_rules! same {
            ($a:expr, $b:expr, $name:literal) => {
                assert!(($a - $b).abs() < 1e-10, "{} changed when an absent locus was added: {} vs {}", $name, $a, $b);
            };
        }

        use crate::calculation::*;
        same!(fst::calculate_fst_for_pair(&d1, &d2, &one, false, 1),
              fst::calculate_fst_for_pair(&d1, &d2, &with_absent, false, 2), "fst");
        same!(gst::calculate_gst_for_pair(&d1, &d2, &one),
              gst::calculate_gst_for_pair(&d1, &d2, &with_absent), "gst");
        same!(chord::calculate_chord_distance_for_pair(&d1, &d2, &one, false, 1),
              chord::calculate_chord_distance_for_pair(&d1, &d2, &with_absent, false, 2), "chord");
        same!(bray_curtis::calculate_bray_curtis_for_pair(&d1, &d2, &one, false, 1),
              bray_curtis::calculate_bray_curtis_for_pair(&d1, &d2, &with_absent, false, 2), "bray-curtis");
        same!(jost_d::calculate_jost_d_for_pair(&d1, &d2, &one, false, 1),
              jost_d::calculate_jost_d_for_pair(&d1, &d2, &with_absent, false, 2), "jost_d");
        same!(reynolds::calculate_reynolds_distance_for_pair(&d1, &d2, &one),
              reynolds::calculate_reynolds_distance_for_pair(&d1, &d2, &with_absent), "reynolds");
    }

    /// A shared monomorphic locus is genetic identity, so it can only pull Nei's D
    /// down, never up.
    #[test]
    fn absent_site_adds_identity_to_nei() {
        let d1 = make_data(&[(100, "A", &[("T", 0.7)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.2)])]);
        let near = crate::calculation::nei::calculate_nei_distance_for_pair(&d1, &d2, &gps(&[100]));
        let far = crate::calculation::nei::calculate_nei_distance_for_pair(&d1, &d2, &gps(&[100, 999]));
        assert!(far < near, "a shared monomorphic locus must reduce Nei's D: {} vs {}", far, near);
        assert!(far > 0.0);
    }

    /// Two identical samples must be at distance 0 no matter how many loci the rest
    /// of the cohort contributes.
    #[test]
    fn identical_samples_zero_despite_foreign_loci() {
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let p = gps(&[100, 200, 300, 400]);
        assert!(crate::calculation::chord::calculate_chord_distance_for_pair(&d1, &d2, &p, false, 4).abs() < 1e-10);
        assert!(crate::calculation::gst::calculate_gst_for_pair(&d1, &d2, &p).abs() < 1e-10);
        assert!(crate::calculation::rogers::calculate_rogers_distance_for_pair(&d1, &d2, &p, 4).abs() < 1e-10);
    }

    #[test]
    fn freqs_multiallelic() {
        let s1 = make_site("A", &[("T", 0.3), ("G", 0.2)]);
        let s2 = make_site("A", &[("T", 0.5)]);
        let psf = get_all_freqs_at_pos(Some(&s1), Some(&s2));
        assert_eq!(psf.all_alleles.len(), 3);
        assert!((psf.freqs1["A"] - 0.5).abs() < 1e-10);
    }

    // ── heterozygosity ──────────────────────────────────────────────────

    #[test]
    fn het_monomorphic() {
        let mut f = AlleleFrequencies::new(); f.insert("A".into(), 1.0);
        assert!((heterozygosity(&f)).abs() < 1e-10);
        assert!((homozygosity(&f) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn het_biallelic() {
        let mut f = AlleleFrequencies::new(); f.insert("A".into(), 0.5); f.insert("T".into(), 0.5);
        assert!((heterozygosity(&f) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn het_plus_hom() {
        let mut f = AlleleFrequencies::new(); f.insert("A".into(), 0.7); f.insert("T".into(), 0.3);
        assert!((heterozygosity(&f) + homozygosity(&f) - 1.0).abs() < 1e-10);
    }

    // ── pooled_heterozygosity ───────────────────────────────────────────

    #[test]
    fn pooled_identical() {
        let s1 = make_site("A", &[("T", 0.3)]);
        let s2 = make_site("A", &[("T", 0.3)]);
        let psf = get_all_freqs_at_pos(Some(&s1), Some(&s2));
        let h_s = (heterozygosity(&psf.freqs1) + heterozygosity(&psf.freqs2)) / 2.0;
        assert!((pooled_heterozygosity(&psf) - h_s).abs() < 1e-10);
    }

    #[test]
    fn pooled_fixed_diff() {
        let s1 = make_site("A", &[("T", 1.0)]);
        let s2 = make_site("A", &[]);
        let psf = get_all_freqs_at_pos(Some(&s1), Some(&s2));
        assert!((pooled_heterozygosity(&psf) - 0.5).abs() < 1e-10);
    }

    // ── FST ─────────────────────────────────────────────────────────────

    #[test]
    fn fst_identical_zero() {
        let d1 = make_data(&[(100, "A", &[("T", 0.3)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.3)])]);
        let fst = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &gps(&[100]), false, 1);
        assert!(fst.abs() < 1e-10);
    }

    #[test]
    fn fst_fixed_one() {
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        let fst = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &gps(&[100]), false, 1);
        assert!((fst - 1.0).abs() < 1e-10);
    }

    #[test]
    fn fst_normalized() {
        let d1 = make_data(&[(100, "A", &[("T", 0.8)]), (200, "C", &[("G", 0.8)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.2)]), (200, "C", &[("G", 0.2)])]);
        let raw = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &gps(&[100, 200]), false, 2);
        let norm = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &gps(&[100, 200]), true, 2);
        assert!((norm - raw / 2.0).abs() < 1e-10);
    }

    #[test]
    fn fst_bounded() {
        for af in &[0.0, 0.1, 0.5, 0.9, 1.0] {
            let d1 = make_data(&[(100, "A", &[("T", *af)])]);
            let d2 = make_data(&[(100, "A", &[("T", 1.0 - af)])]);
            let fst = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &gps(&[100]), false, 1);
            assert!(fst >= -1e-10 && fst <= 1.0 + 1e-10);
        }
    }

    // ── GST ─────────────────────────────────────────────────────────────

    #[test]
    fn gst_identical_zero() {
        let d1 = make_data(&[(100, "A", &[("T", 0.5)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.5)])]);
        assert!(crate::calculation::gst::calculate_gst_for_pair(&d1, &d2, &gps(&[100])).abs() < 1e-10);
    }

    #[test]
    fn gst_fixed_one() {
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        assert!((crate::calculation::gst::calculate_gst_for_pair(&d1, &d2, &gps(&[100])) - 1.0).abs() < 1e-10);
    }

    // ── Jost's D ────────────────────────────────────────────────────────

    #[test]
    fn jost_d_identical_zero() {
        let d1 = make_data(&[(100, "A", &[("T", 0.3)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.3)])]);
        assert!(crate::calculation::jost_d::calculate_jost_d_for_pair(&d1, &d2, &gps(&[100]), false, 1).abs() < 1e-10);
    }

    #[test]
    fn jost_d_fixed_one() {
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        assert!((crate::calculation::jost_d::calculate_jost_d_for_pair(&d1, &d2, &gps(&[100]), false, 1) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn jost_d_non_negative() {
        for af1 in &[0.0, 0.3, 0.5, 0.7, 1.0] {
            for af2 in &[0.0, 0.3, 0.5, 0.7, 1.0] {
                let d1 = make_data(&[(100, "A", &[("T", *af1)])]);
                let d2 = make_data(&[(100, "A", &[("T", *af2)])]);
                let d = crate::calculation::jost_d::calculate_jost_d_for_pair(&d1, &d2, &gps(&[100]), false, 1);
                assert!(d >= -1e-10);
            }
        }
    }

    #[test]
    fn jost_d_symmetric() {
        let d1 = make_data(&[(100, "A", &[("T", 0.7)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.2)])]);
        let a = crate::calculation::jost_d::calculate_jost_d_for_pair(&d1, &d2, &gps(&[100]), false, 1);
        let b = crate::calculation::jost_d::calculate_jost_d_for_pair(&d2, &d1, &gps(&[100]), false, 1);
        assert!((a - b).abs() < 1e-10);
    }

    // ── Nei ─────────────────────────────────────────────────────────────

    #[test]
    fn nei_identical_zero() {
        let d1 = make_data(&[(100, "A", &[("T", 0.4)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.4)])]);
        assert!(crate::calculation::nei::calculate_nei_distance_for_pair(&d1, &d2, &gps(&[100])).abs() < 1e-10);
    }

    #[test]
    fn nei_fixed_inf() {
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        assert!(crate::calculation::nei::calculate_nei_distance_for_pair(&d1, &d2, &gps(&[100])).is_infinite());
    }

    #[test]
    fn nei_symmetric() {
        let d1 = make_data(&[(100, "A", &[("T", 0.8)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.3)])]);
        let a = crate::calculation::nei::calculate_nei_distance_for_pair(&d1, &d2, &gps(&[100]));
        let b = crate::calculation::nei::calculate_nei_distance_for_pair(&d2, &d1, &gps(&[100]));
        assert!((a - b).abs() < 1e-10);
    }

    // ── Chord ───────────────────────────────────────────────────────────

    #[test]
    fn chord_identical_zero() {
        let d1 = make_data(&[(100, "A", &[("T", 0.5)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.5)])]);
        assert!(crate::calculation::chord::calculate_chord_distance_for_pair(&d1, &d2, &gps(&[100]), false, 1).abs() < 1e-10);
    }

    #[test]
    fn chord_fixed() {
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        let d = crate::calculation::chord::calculate_chord_distance_for_pair(&d1, &d2, &gps(&[100]), false, 1);
        assert!((d - 2.0_f64.sqrt()).abs() < 1e-10);
    }

    // ── Bray-Curtis ─────────────────────────────────────────────────────

    #[test]
    fn bc_identical_zero() {
        let d1 = make_data(&[(100, "A", &[("T", 0.3)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.3)])]);
        assert!(crate::calculation::bray_curtis::calculate_bray_curtis_for_pair(&d1, &d2, &gps(&[100]), false, 1).abs() < 1e-10);
    }

    #[test]
    fn bc_fixed() {
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        assert!((crate::calculation::bray_curtis::calculate_bray_curtis_for_pair(&d1, &d2, &gps(&[100]), false, 1) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn bc_bounded() {
        for af in &[0.0, 0.2, 0.5, 0.8, 1.0] {
            let d1 = make_data(&[(100, "A", &[("T", *af)])]);
            let d2 = make_data(&[(100, "A", &[("T", 1.0 - af)])]);
            let d = crate::calculation::bray_curtis::calculate_bray_curtis_for_pair(&d1, &d2, &gps(&[100]), false, 1);
            assert!(d >= -1e-10 && d <= 1.0 + 1e-10);
        }
    }

    // ── Reynolds ────────────────────────────────────────────────────────

    #[test]
    fn reynolds_identical_zero() {
        let d1 = make_data(&[(100, "A", &[("T", 0.5)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.5)])]);
        assert!(crate::calculation::reynolds::calculate_reynolds_distance_for_pair(&d1, &d2, &gps(&[100])).abs() < 1e-10);
    }

    #[test]
    fn reynolds_fixed_inf() {
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        assert!(crate::calculation::reynolds::calculate_reynolds_distance_for_pair(&d1, &d2, &gps(&[100])).is_infinite());
    }

    /// At intermediate frequencies the coancestry estimator and Nei's GST diverge, so
    /// this is the case that tells them apart. Substituting GST would give 0.4054651.
    #[test]
    fn reynolds_uses_coancestry_not_gst() {
        // p = {T:0.5, A:0.5} vs q = {T:1.0, A:0.0}
        // theta = ((0.5)^2 + (0.5)^2)/2 / (1 - 0.5) = 0.25/0.5 = 0.5
        let d1 = make_data(&[(100, "A", &[("T", 0.5)])]);
        let d2 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d = crate::calculation::reynolds::calculate_reynolds_distance_for_pair(&d1, &d2, &gps(&[100]));
        assert!((d - (2.0_f64).ln()).abs() < 1e-10, "expected -ln(0.5), got {}", d);
    }

    #[test]
    fn reynolds_symmetric() {
        let d1 = make_data(&[(100, "A", &[("T", 0.7)]), (200, "C", &[("G", 0.1)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.2)]), (200, "C", &[("G", 0.9)])]);
        let p = gps(&[100, 200]);
        let a = crate::calculation::reynolds::calculate_reynolds_distance_for_pair(&d1, &d2, &p);
        let b = crate::calculation::reynolds::calculate_reynolds_distance_for_pair(&d2, &d1, &p);
        assert!((a - b).abs() < 1e-10);
        assert!(a > 0.0);
    }

    // ── Rogers ──────────────────────────────────────────────────────────

    #[test]
    fn rogers_identical_zero() {
        let d1 = make_data(&[(100, "A", &[("T", 0.5)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.5)])]);
        assert!(crate::calculation::rogers::calculate_rogers_distance_for_pair(&d1, &d2, &gps(&[100]), 1).abs() < 1e-10);
    }

    #[test]
    fn rogers_fixed() {
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        assert!((crate::calculation::rogers::calculate_rogers_distance_for_pair(&d1, &d2, &gps(&[100]), 1) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn rogers_bounded() {
        for af in &[0.0, 0.2, 0.5, 0.8, 1.0] {
            let d1 = make_data(&[(100, "A", &[("T", *af)])]);
            let d2 = make_data(&[(100, "A", &[("T", 1.0 - af)])]);
            let d = crate::calculation::rogers::calculate_rogers_distance_for_pair(&d1, &d2, &gps(&[100]), 1);
            assert!(d >= -1e-10 && d <= 1.0 + 1e-10);
        }
    }

    // ── Utility ─────────────────────────────────────────────────────────

    #[test]
    fn pairs_count() {
        let s: Vec<String> = (0..5).map(|i| format!("s{}", i)).collect();
        assert_eq!(crate::calculation::generate_sample_pairs(&s).len(), 10);
    }

    #[test]
    fn pairs_empty() { assert!(crate::calculation::generate_sample_pairs(&[]).is_empty()); }

    #[test]
    fn pairs_single() { assert!(crate::calculation::generate_sample_pairs(&["s1".into()]).is_empty()); }

    #[test]
    fn matrix_symmetric() {
        let r = vec![((0, 1), 0.5), ((0, 2), 0.3), ((1, 2), 0.7)];
        let m = crate::calculation::create_distance_matrix(&r, 3);

        // m[0][1] == m[1][0] holds even if both cells hold the wrong number, so check
        // the values that went in actually came back out.
        assert_eq!(m[0][1], 0.5);
        assert_eq!(m[1][0], 0.5);
        assert_eq!(m[0][2], 0.3);
        assert_eq!(m[2][0], 0.3);
        assert_eq!(m[1][2], 0.7);
        assert_eq!(m[2][1], 0.7);
        for (i, row) in m.iter().enumerate() {
            assert_eq!(row[i], 0.0, "diagonal at {}", i);
        }
    }

    #[test]
    fn quote_field_protects_separators() {
        use crate::io::csv::write_distance_matrix;
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.csv");
        let m = vec![vec![0.0, 0.5], vec![0.5, 0.0]];
        let samples = vec!["with,comma".to_string(), "he said \"hi\"".to_string()];
        write_distance_matrix(out.to_str().unwrap(), &m, &samples, false).unwrap();

        let c = std::fs::read_to_string(&out).unwrap();
        assert!(c.contains("\"with,comma\""), "got: {}", c);
        // Embedded quotes are doubled, per RFC 4180
        assert!(c.contains("\"he said \"\"hi\"\"\""), "got: {}", c);
    }

    #[test]
    fn non_finite_cells_written_as_na() {
        use crate::io::csv::write_distance_matrix;
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.csv");
        let m = vec![
            vec![0.0, f64::INFINITY],
            vec![f64::INFINITY, 0.0],
        ];
        write_distance_matrix(out.to_str().unwrap(), &m, &["A".into(), "B".into()], false).unwrap();

        let c = std::fs::read_to_string(&out).unwrap();
        assert!(c.contains("A,0.0000000000,NA"), "got: {}", c);
        assert!(!c.to_lowercase().contains("inf"), "got: {}", c);
    }

    #[test]
    fn fst_multi_locus() {
        let d1 = make_data(&[(100, "A", &[("T", 0.9)]), (200, "C", &[("G", 0.1)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.1)]), (200, "C", &[("G", 0.9)])]);
        let fst = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &gps(&[100, 200]), false, 2);
        let single = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &gps(&[100]), false, 1);
        assert!((fst - 2.0 * single).abs() < 1e-10);
    }

    // ── Multi-chrom correctness ─────────────────────────────────────────

    #[test]
    fn different_chroms_are_distinct_sites() {
        let pos_a = GenomicPos { chrom: "chr1".into(), pos: 100 };
        let pos_b = GenomicPos { chrom: "chr2".into(), pos: 100 };
        assert_ne!(pos_a, pos_b);

        // Same position number on different chroms should be independent
        let mut d1 = PositionalData::new();
        d1.insert(pos_a.clone(), make_site("A", &[("T", 0.9)]));
        d1.insert(pos_b.clone(), make_site("C", &[("G", 0.1)]));

        let mut d2 = PositionalData::new();
        d2.insert(pos_a.clone(), make_site("A", &[("T", 0.1)]));
        d2.insert(pos_b.clone(), make_site("C", &[("G", 0.9)]));

        let positions = vec![pos_a, pos_b];
        let fst = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &positions, false, 2);
        assert!(fst > 0.0, "Different chroms should contribute independently");
    }

    // ── FASTA ───────────────────────────────────────────────────────────

    #[test]
    fn fasta_multi_contig() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("ref.fa");
        std::fs::write(&p, ">chr1\nACGT\n>chr2\nTTTT\n").unwrap();
        let g = crate::io::fasta::read_reference(p.to_str().unwrap()).unwrap();
        assert_eq!(g.len(), 2);
        assert_eq!(g["chr1"], b"ACGT");
        assert_eq!(g["chr2"], b"TTTT");
    }

    #[test]
    fn fasta_deterministic() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("ref.fa");
        std::fs::write(&p, ">z_last\nGGGG\n>a_first\nACGT\n").unwrap();
        let g = crate::io::fasta::read_reference(p.to_str().unwrap()).unwrap();
        let names: Vec<&str> = g.keys().map(|s| s.as_str()).collect();
        assert_eq!(names, vec!["a_first", "z_last"]);
    }

    #[test]
    fn fasta_softmasked_uppercase() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("ref.fa");
        std::fs::write(&p, ">chr1\nacgt\nNNnn\n").unwrap();
        let g = crate::io::fasta::read_reference(p.to_str().unwrap()).unwrap();
        assert_eq!(g["chr1"], b"ACGTNNNN");
    }

    #[test]
    fn fasta_empty() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("empty.fa");
        std::fs::write(&p, "").unwrap();
        assert!(crate::io::fasta::read_reference(p.to_str().unwrap()).is_err());
    }

    fn read_fa(body: &str) -> Result<crate::io::fasta::ReferenceGenome, std::io::Error> {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("ref.fa");
        std::fs::write(&p, body).unwrap();
        crate::io::fasta::read_reference(p.to_str().unwrap())
    }

    #[test]
    fn fasta_nameless_record_rejected() {
        // Used to leave current_name empty, so the nameless record's bases were
        // prepended to the next contig and shifted every coordinate after it.
        assert!(read_fa(">\nAAAA\n>chr1\nCCCC\n").is_err());
    }

    #[test]
    fn fasta_sequence_before_header_rejected() {
        assert!(read_fa("AAAA\n>chr1\nCCCC\n").is_err());
    }

    #[test]
    fn fasta_duplicate_contig_rejected() {
        assert!(read_fa(">chr1\nAAAA\n>chr1\nCCCC\n").is_err());
    }

    #[test]
    fn fasta_records_stay_separate() {
        let g = read_fa(">chr1 description here\nAA\nAA\n>chr2\nCCCC\n").unwrap();
        assert_eq!(g["chr1"], b"AAAA");
        assert_eq!(g["chr2"], b"CCCC");
    }

    // ── Frequency renormalisation ───────────────────────────────────────

    #[test]
    fn saturated_site_is_rescaled() {
        // Split multi-allelic records: two independently estimated 0.9 frequencies.
        let mut sv = crate::types::SampleVariants::new();
        sv.insert("s1".to_string(), make_data(&[(100, "A", &[("T", 0.9), ("G", 0.9)])]));
        assert_eq!(crate::io::renormalise_saturated_sites(&mut sv), 1);

        let freqs = &sv["s1"][&gp(100)].freqs;
        let sum: f64 = freqs.values().sum();
        assert!((sum - 1.0).abs() < 1e-10, "sum is {}", sum);
        assert!((freqs["T"] - 0.5).abs() < 1e-10);
        // Heterozygosity must be back in [0,1]
        let h = heterozygosity(freqs);
        assert!((0.0..=1.0).contains(&h), "heterozygosity is {}", h);
    }

    #[test]
    fn well_formed_sites_are_left_alone() {
        let mut sv = crate::types::SampleVariants::new();
        sv.insert("s1".to_string(), make_data(&[(100, "A", &[("T", 0.3)]), (200, "C", &[("G", 1.0)])]));
        assert_eq!(crate::io::renormalise_saturated_sites(&mut sv), 0);
        assert!((sv["s1"][&gp(100)].freqs["T"] - 0.3).abs() < 1e-10);
    }

    #[test]
    fn rounding_noise_is_not_rescaled() {
        let mut sv = crate::types::SampleVariants::new();
        sv.insert("s1".to_string(), make_data(&[(100, "A", &[("T", 1.0 + 1e-12)])]));
        assert_eq!(crate::io::renormalise_saturated_sites(&mut sv), 0);
    }

    // ── All symmetric ───────────────────────────────────────────────────

    #[test]
    fn all_metrics_symmetric() {
        let d1 = make_data(&[(100, "A", &[("T", 0.7)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.2)])]);
        let p = gps(&[100]);

        let f12 = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &p, false, 1);
        let f21 = crate::calculation::fst::calculate_fst_for_pair(&d2, &d1, &p, false, 1);
        assert!((f12 - f21).abs() < 1e-10);

        let b12 = crate::calculation::bray_curtis::calculate_bray_curtis_for_pair(&d1, &d2, &p, false, 1);
        let b21 = crate::calculation::bray_curtis::calculate_bray_curtis_for_pair(&d2, &d1, &p, false, 1);
        assert!((b12 - b21).abs() < 1e-10);

        let c12 = crate::calculation::chord::calculate_chord_distance_for_pair(&d1, &d2, &p, false, 1);
        let c21 = crate::calculation::chord::calculate_chord_distance_for_pair(&d2, &d1, &p, false, 1);
        assert!((c12 - c21).abs() < 1e-10);

        let r12 = crate::calculation::rogers::calculate_rogers_distance_for_pair(&d1, &d2, &p, 1);
        let r21 = crate::calculation::rogers::calculate_rogers_distance_for_pair(&d2, &d1, &p, 1);
        assert!((r12 - r21).abs() < 1e-10);
    }

    // ── Output format ───────────────────────────────────────────────────

    #[test]
    fn write_csv() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.csv");
        let m = vec![vec![0.0, 0.5], vec![0.5, 0.0]];
        crate::io::csv::write_distance_matrix(out.to_str().unwrap(), &m, &["A".into(), "B".into()], false).unwrap();
        let c = std::fs::read_to_string(&out).unwrap();
        assert!(c.starts_with("sample,"));
        assert!(c.contains("A,0.0000000000,0.5000000000"));
    }

    #[test]
    fn write_tsv() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("out.tsv");
        let m = vec![vec![0.0, 0.5], vec![0.5, 0.0]];
        crate::io::csv::write_distance_matrix(out.to_str().unwrap(), &m, &["A".into(), "B".into()], true).unwrap();
        let c = std::fs::read_to_string(&out).unwrap();
        assert!(c.starts_with("sample\t"));
    }

    // ── Ref imputation does not overwrite explicit freq ──────────────

    #[test]
    fn ref_imputation_no_overwrite() {
        // If a sample explicitly has the ref allele as an alt with freq 0.7,
        // imputation must NOT overwrite it with 1-sum.
        let s1 = make_site("A", &[("A", 0.7), ("T", 0.3)]);
        let s2 = make_site("A", &[("T", 0.5)]);
        let psf = get_all_freqs_at_pos(Some(&s1), Some(&s2));
        // s1 already has "A" in freqs → should keep 0.7, not overwrite with 1-1.0=0.0
        assert!((psf.freqs1["A"] - 0.7).abs() < 1e-10);
        // s2 does NOT have "A" in freqs → should impute 1-0.5=0.5
        assert!((psf.freqs2["A"] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn ref_imputation_normal_case() {
        // Normal VCF case: only alt alleles in freqs, ref should be imputed
        let s1 = make_site("A", &[("T", 0.3)]);
        let s2 = make_site("A", &[("T", 0.6)]);
        let psf = get_all_freqs_at_pos(Some(&s1), Some(&s2));
        assert!((psf.freqs1["A"] - 0.7).abs() < 1e-10);
        assert!((psf.freqs2["A"] - 0.4).abs() < 1e-10);
    }
}
