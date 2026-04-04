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
        let s1 = SiteData { reference_allele: String::new(), freqs: { let mut f = AlleleFrequencies::new(); f.insert("T".into(), 0.5); f } };
        let psf = get_all_freqs_at_pos(Some(&s1), None);
        assert_eq!(psf.freqs1.len(), 1);
    }

    #[test]
    fn freqs_both_none() {
        let psf = get_all_freqs_at_pos(None, None);
        assert!(psf.all_alleles.is_empty());
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
        assert_eq!(m[0][1], m[1][0]);
        assert_eq!(m[0][0], 0.0);
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
        assert_eq!(crate::io::fasta::get_base_at(&g, 1), Some(b'A'));
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
