#[cfg(test)]
mod tests {
    use crate::calculation::common::*;
    use crate::types::{AlleleFrequencies, PositionalData, SiteData};
    use std::collections::{HashMap, HashSet};
    extern crate tempfile;

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
        let site = make_site("A", &[("T", 0.7), ("G", 0.4)]);
        let psf = get_all_freqs_at_pos(Some(&site), None);
        assert!(psf.freqs1["A"] >= 0.0);
    }

    #[test]
    fn freqs_no_ref_allele() {
        let s1 = SiteData {
            reference_allele: String::new(),
            freqs: {
                let mut f = AlleleFrequencies::new();
                f.insert("T".into(), 0.5);
                f
            },
        };
        let psf = get_all_freqs_at_pos(Some(&s1), None);
        assert!(psf.all_alleles.contains("T"));
        // No ref allele → freqs should only have what's explicit
        assert_eq!(psf.freqs1.len(), 1);
    }

    #[test]
    fn freqs_both_none() {
        let psf = get_all_freqs_at_pos(None, None);
        assert!(psf.all_alleles.is_empty());
        assert!(psf.freqs1.is_empty());
        assert!(psf.freqs2.is_empty());
    }

    #[test]
    fn freqs_multiallelic() {
        let s1 = make_site("A", &[("T", 0.3), ("G", 0.2)]);
        let s2 = make_site("A", &[("T", 0.5)]);
        let psf = get_all_freqs_at_pos(Some(&s1), Some(&s2));
        assert_eq!(psf.all_alleles.len(), 3); // A, T, G
        assert!((psf.freqs1["A"] - 0.5).abs() < 1e-10);
        assert!((psf.freqs2.get("G").copied().unwrap_or(0.0)).abs() < 1e-10); // not present in s2
    }

    // ── heterozygosity / homozygosity ───────────────────────────────────

    #[test]
    fn het_monomorphic() {
        let mut f = AlleleFrequencies::new();
        f.insert("A".into(), 1.0);
        assert!((heterozygosity(&f) - 0.0).abs() < 1e-10);
        assert!((homozygosity(&f) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn het_biallelic_equal() {
        let mut f = AlleleFrequencies::new();
        f.insert("A".into(), 0.5);
        f.insert("T".into(), 0.5);
        assert!((heterozygosity(&f) - 0.5).abs() < 1e-10);
        assert!((homozygosity(&f) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn het_triallelic() {
        let mut f = AlleleFrequencies::new();
        f.insert("A".into(), 1.0 / 3.0);
        f.insert("T".into(), 1.0 / 3.0);
        f.insert("G".into(), 1.0 / 3.0);
        // H = 1 - 3*(1/3)^2 = 1 - 1/3 = 2/3
        assert!((heterozygosity(&f) - 2.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn het_plus_hom_equals_one() {
        let mut f = AlleleFrequencies::new();
        f.insert("A".into(), 0.7);
        f.insert("T".into(), 0.3);
        let h = heterozygosity(&f);
        let j = homozygosity(&f);
        assert!((h + j - 1.0).abs() < 1e-10);
    }

    // ── pooled_heterozygosity ───────────────────────────────────────────

    #[test]
    fn pooled_het_identical_pops() {
        let s1 = make_site("A", &[("T", 0.3)]);
        let s2 = make_site("A", &[("T", 0.3)]);
        let psf = get_all_freqs_at_pos(Some(&s1), Some(&s2));
        let h_s = (heterozygosity(&psf.freqs1) + heterozygosity(&psf.freqs2)) / 2.0;
        let h_t = pooled_heterozygosity(&psf);
        assert!((h_t - h_s).abs() < 1e-10, "Identical pops: Ht should equal Hs");
    }

    #[test]
    fn pooled_het_fixed_diff() {
        let s1 = make_site("A", &[("T", 1.0)]); // 100% T
        let s2 = make_site("A", &[]); // 100% A (ref)
        let psf = get_all_freqs_at_pos(Some(&s1), Some(&s2));
        let h_t = pooled_heterozygosity(&psf);
        assert!((h_t - 0.5).abs() < 1e-10, "Fixed diff: Ht should be 0.5");
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
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        let pos = positions(&[100]);
        let fst = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &pos, false, 1);
        assert!((fst - 1.0).abs() < 1e-10, "Fixed diff → FST=1, got {}", fst);
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

    #[test]
    fn fst_bounded_0_1() {
        for af in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
            let d1 = make_data(&[(100, "A", &[("T", *af)])]);
            let d2 = make_data(&[(100, "A", &[("T", 1.0 - af)])]);
            let pos = positions(&[100]);
            let fst = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &pos, false, 1);
            assert!(fst >= -1e-10 && fst <= 1.0 + 1e-10, "FST should be in [0,1], got {} for af={}", fst, af);
        }
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

    #[test]
    fn gst_fixed_diff_is_one() {
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        let pos = positions(&[100]);
        let gst = crate::calculation::gst::calculate_gst_for_pair(&d1, &d2, &pos);
        assert!((gst - 1.0).abs() < 1e-10);
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
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        let pos = positions(&[100]);
        let d = crate::calculation::jost_d::calculate_jost_d_for_pair(&d1, &d2, &pos, false, 1);
        // Hs=0, Ht=0.5 → D = 2*(0.5-0)/(1-0) = 1.0
        assert!((d - 1.0).abs() < 1e-10, "Fixed diff → D=1, got {}", d);
    }

    #[test]
    fn jost_d_is_always_non_negative() {
        for af1 in &[0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0] {
            for af2 in &[0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0] {
                let d1 = make_data(&[(100, "A", &[("T", *af1)])]);
                let d2 = make_data(&[(100, "A", &[("T", *af2)])]);
                let pos = positions(&[100]);
                let d = crate::calculation::jost_d::calculate_jost_d_for_pair(&d1, &d2, &pos, false, 1);
                assert!(d >= -1e-10, "Jost D >= 0, got {} for af1={} af2={}", d, af1, af2);
            }
        }
    }

    #[test]
    fn jost_d_symmetrical() {
        let d1 = make_data(&[(100, "A", &[("T", 0.7)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.2)])]);
        let pos = positions(&[100]);
        let d_12 = crate::calculation::jost_d::calculate_jost_d_for_pair(&d1, &d2, &pos, false, 1);
        let d_21 = crate::calculation::jost_d::calculate_jost_d_for_pair(&d2, &d1, &pos, false, 1);
        assert!((d_12 - d_21).abs() < 1e-10, "Jost D should be symmetric");
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

    #[test]
    fn nei_symmetric() {
        let d1 = make_data(&[(100, "A", &[("T", 0.8)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.3)])]);
        let pos = positions(&[100]);
        let d_12 = crate::calculation::nei::calculate_nei_distance_for_pair(&d1, &d2, &pos);
        let d_21 = crate::calculation::nei::calculate_nei_distance_for_pair(&d2, &d1, &pos);
        assert!((d_12 - d_21).abs() < 1e-10);
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
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn bray_curtis_bounded() {
        for af in &[0.0, 0.2, 0.5, 0.8, 1.0] {
            let d1 = make_data(&[(100, "A", &[("T", *af)])]);
            let d2 = make_data(&[(100, "A", &[("T", 1.0 - af)])]);
            let pos = positions(&[100]);
            let d = crate::calculation::bray_curtis::calculate_bray_curtis_for_pair(&d1, &d2, &pos, false, 1);
            assert!(d >= -1e-10 && d <= 1.0 + 1e-10, "BC in [0,1], got {}", d);
        }
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

    #[test]
    fn reynolds_fixed_diff_is_inf() {
        let d1 = make_data(&[(100, "A", &[("T", 1.0)])]);
        let d2 = make_data(&[(100, "A", &[])]);
        let pos = positions(&[100]);
        let d = crate::calculation::reynolds::calculate_reynolds_distance_for_pair(&d1, &d2, &pos);
        assert!(d.is_infinite(), "Fixed diff → inf, got {}", d);
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
        assert!((d - 1.0).abs() < 1e-10);
    }

    #[test]
    fn rogers_bounded_0_1() {
        for af in &[0.0, 0.2, 0.5, 0.8, 1.0] {
            let d1 = make_data(&[(100, "A", &[("T", *af)])]);
            let d2 = make_data(&[(100, "A", &[("T", 1.0 - af)])]);
            let pos = positions(&[100]);
            let d = crate::calculation::rogers::calculate_rogers_distance_for_pair(&d1, &d2, &pos, 1);
            assert!(d >= -1e-10 && d <= 1.0 + 1e-10, "Rogers in [0,1], got {}", d);
        }
    }

    // ── generate_sample_pairs ───────────────────────────────────────────

    #[test]
    fn pairs_count() {
        let samples: Vec<String> = (0..5).map(|i| format!("s{}", i)).collect();
        let pairs = crate::calculation::generate_sample_pairs(&samples);
        assert_eq!(pairs.len(), 10);
    }

    #[test]
    fn pairs_empty() {
        let samples: Vec<String> = vec![];
        let pairs = crate::calculation::generate_sample_pairs(&samples);
        assert!(pairs.is_empty());
    }

    #[test]
    fn pairs_single_sample() {
        let samples = vec!["s1".to_string()];
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

    #[test]
    fn matrix_diagonal_zero() {
        let results = vec![((0, 1), 0.5)];
        let m = crate::calculation::create_distance_matrix(&results, 2);
        assert_eq!(m[0][0], 0.0);
        assert_eq!(m[1][1], 0.0);
    }

    // ── multi-locus consistency ──────────────────────────────────────────

    #[test]
    fn fst_multi_locus_sum() {
        let d1 = make_data(&[(100, "A", &[("T", 0.9)]), (200, "C", &[("G", 0.1)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.1)]), (200, "C", &[("G", 0.9)])]);
        let pos = positions(&[100, 200]);
        let fst = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &pos, false, 2);
        let single_pos = positions(&[100]);
        let fst_single = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &single_pos, false, 1);
        assert!((fst - 2.0 * fst_single).abs() < 1e-10);
    }

    // ── FASTA reader ────────────────────────────────────────────────────

    #[test]
    fn fasta_multi_contig() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("ref.fa");
        std::fs::write(&fasta, ">chr1\nACGT\n>chr2\nTTTT\n").unwrap();
        let genome = crate::io::fasta::read_reference(fasta.to_str().unwrap()).unwrap();
        assert_eq!(genome.len(), 2);
        assert_eq!(genome["chr1"], b"ACGT");
        assert_eq!(genome["chr2"], b"TTTT");
    }

    #[test]
    fn fasta_multiline_seq() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("ref.fa");
        std::fs::write(&fasta, ">contig1\nACGT\nAAAA\n").unwrap();
        let genome = crate::io::fasta::read_reference(fasta.to_str().unwrap()).unwrap();
        assert_eq!(genome["contig1"], b"ACGTAAAA");
    }

    #[test]
    fn fasta_get_base() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("ref.fa");
        std::fs::write(&fasta, ">chr1\nACGT\n").unwrap();
        let genome = crate::io::fasta::read_reference(fasta.to_str().unwrap()).unwrap();
        assert_eq!(crate::io::fasta::get_base_at(&genome, 1), Some(b'A'));
        assert_eq!(crate::io::fasta::get_base_at(&genome, 4), Some(b'T'));
        assert_eq!(crate::io::fasta::get_base_at(&genome, 5), None);
    }

    #[test]
    fn fasta_empty_file() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("empty.fa");
        std::fs::write(&fasta, "").unwrap();
        let result = crate::io::fasta::read_reference(fasta.to_str().unwrap());
        assert!(result.is_err());
    }

    // ── All metrics symmetric ───────────────────────────────────────────

    #[test]
    fn all_metrics_symmetric() {
        let d1 = make_data(&[(100, "A", &[("T", 0.7)])]);
        let d2 = make_data(&[(100, "A", &[("T", 0.2)])]);
        let pos = positions(&[100]);

        let fst_12 = crate::calculation::fst::calculate_fst_for_pair(&d1, &d2, &pos, false, 1);
        let fst_21 = crate::calculation::fst::calculate_fst_for_pair(&d2, &d1, &pos, false, 1);
        assert!((fst_12 - fst_21).abs() < 1e-10, "FST not symmetric");

        let bc_12 = crate::calculation::bray_curtis::calculate_bray_curtis_for_pair(&d1, &d2, &pos, false, 1);
        let bc_21 = crate::calculation::bray_curtis::calculate_bray_curtis_for_pair(&d2, &d1, &pos, false, 1);
        assert!((bc_12 - bc_21).abs() < 1e-10, "BC not symmetric");

        let ch_12 = crate::calculation::chord::calculate_chord_distance_for_pair(&d1, &d2, &pos, false, 1);
        let ch_21 = crate::calculation::chord::calculate_chord_distance_for_pair(&d2, &d1, &pos, false, 1);
        assert!((ch_12 - ch_21).abs() < 1e-10, "Chord not symmetric");

        let ro_12 = crate::calculation::rogers::calculate_rogers_distance_for_pair(&d1, &d2, &pos, 1);
        let ro_21 = crate::calculation::rogers::calculate_rogers_distance_for_pair(&d2, &d1, &pos, 1);
        assert!((ro_12 - ro_21).abs() < 1e-10, "Rogers not symmetric");
    }
}
