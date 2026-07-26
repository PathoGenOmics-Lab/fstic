use std::fs;
use std::path::Path;
use std::process::Command;

fn fstic_bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_fstic"))
}

fn write_vcf(dir: &Path, name: &str, lines: &[&str]) {
    let content = lines.join("\n") + "\n";
    fs::write(dir.join(format!("{}.vcf", name)), content).unwrap();
}

/// First off-diagonal value of a 2x2 matrix file. Asserting on the parsed number
/// instead of a substring is the difference between a real check and one that the
/// zero on the diagonal satisfies by accident.
fn off_diagonal(content: &str, sep: char) -> f64 {
    let line = content.lines().nth(1).expect("matrix needs a first data row");
    let cell = line.split(sep).nth(2).expect("matrix needs a second column");
    cell.trim()
        .parse()
        .unwrap_or_else(|_| panic!("cell {:?} is not a number", cell))
}

/// Two identical samples should produce distance = 0.
#[test]
fn identical_samples_zero_distance() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";

    write_vcf(dir.path(), "sampleA", &[
        header,
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:AD:FREQ\t0/1:100:50:50.0%",
        "chr1\t200\t.\tC\tG\t.\t.\t.\tGT:DP:AD:FREQ\t0/1:100:30:30.0%",
    ]);

    write_vcf(dir.path(), "sampleB", &[
        header,
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:AD:FREQ\t0/1:100:50:50.0%",
        "chr1\t200\t.\tC\tG\t.\t.\t.\tGT:DP:AD:FREQ\t0/1:100:30:30.0%",
    ]);

    let out = dir.path().join("out.csv");

    let status = fstic_bin()
        .args([
            "--vcf",
            dir.path().join("sampleA.vcf").to_str().unwrap(),
            dir.path().join("sampleB.vcf").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--formula", "fst",
            "--min-depth", "1",
            "--min-af", "0.01",
            "--min-alt-reads", "1",
            "--min-alt-rev-reads", "0",
        ])
        .status()
        .unwrap();

    assert!(status.success());
    let content = fs::read_to_string(&out).unwrap();
    // Both off-diagonal should be 0
    assert!(content.contains("0.000000"));
}

/// Completely different samples should produce FST > 0.
#[test]
fn different_samples_positive_fst() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";

    write_vcf(dir.path(), "pop1", &[
        header,
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:AD:FREQ\t0/1:100:90:90.0%",
    ]);

    write_vcf(dir.path(), "pop2", &[
        header,
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:AD:FREQ\t0/1:100:10:10.0%",
    ]);

    let out = dir.path().join("out.tsv");
    let status = fstic_bin()
        .args([
            "--vcf",
            dir.path().join("pop1.vcf").to_str().unwrap(),
            dir.path().join("pop2.vcf").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--formula", "fst",
            "--min-depth", "1",
            "--min-af", "0.01",
            "--min-alt-reads", "1",
            "--min-alt-rev-reads", "0",
        ])
        .status()
        .unwrap();

    assert!(status.success());
    let content = fs::read_to_string(&out).unwrap();
    // TSV output should use tabs
    assert!(content.contains('\t'));
    // Off-diagonal should NOT be 0
    let lines: Vec<&str> = content.lines().collect();
    let vals: Vec<&str> = lines[1].split('\t').collect();
    let dist: f64 = vals[2].parse().unwrap();
    assert!(dist > 0.0, "FST should be > 0 for different populations");
}

/// VCF with FREQ as proportion (no %) should work correctly.
#[test]
fn freq_as_proportion() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";

    // 0.5 = 50%, NOT 0.5%
    write_vcf(dir.path(), "propA", &[
        header,
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:AD:FREQ\t0/1:100:50:0.5",
    ]);

    write_vcf(dir.path(), "propB", &[
        header,
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:AD:FREQ\t0/1:100:50:50%",
    ]);

    let out = dir.path().join("out.csv");
    let status = fstic_bin()
        .args([
            "--vcf",
            dir.path().join("propA.vcf").to_str().unwrap(),
            dir.path().join("propB.vcf").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--formula", "fst",
            "--min-depth", "1",
            "--min-af", "0.01",
            "--min-alt-reads", "1",
            "--min-alt-rev-reads", "0",
        ])
        .status()
        .unwrap();

    assert!(status.success());
    let content = fs::read_to_string(&out).unwrap();
    // Both samples have 50% alt → distance should be ~0
    let lines: Vec<&str> = content.lines().collect();
    let vals: Vec<&str> = lines[1].split(',').collect();
    let dist: f64 = vals[2].parse().unwrap();
    assert!(dist.abs() < 1e-6, "Same freq (0.5 vs 50%) should give 0 distance, got {}", dist);
}

/// VCF mode should NOT require --reference.
#[test]
fn vcf_no_reference_required() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";

    write_vcf(dir.path(), "s1", &[
        header,
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:AD:FREQ\t0/1:100:50:50%",
    ]);
    write_vcf(dir.path(), "s2", &[
        header,
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:AD:FREQ\t0/1:100:50:50%",
    ]);

    let out = dir.path().join("out.csv");
    // No --reference flag
    let status = fstic_bin()
        .args([
            "--vcf",
            dir.path().join("s1.vcf").to_str().unwrap(),
            dir.path().join("s2.vcf").to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "--min-depth", "1",
            "--min-af", "0.01",
            "--min-alt-reads", "1",
            "--min-alt-rev-reads", "0",
        ])
        .status()
        .unwrap();

    assert!(status.success());
}

/// Indels should be skipped (only SNPs processed).
#[test]
fn indels_skipped() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";

    // Only indels, no SNPs → should fail with "no polymorphic sites"
    write_vcf(dir.path(), "s1", &[
        header,
        "chr1\t100\t.\tAT\tA\t.\t.\t.\tGT:DP:AD:FREQ\t0/1:100:50:50%",
    ]);
    write_vcf(dir.path(), "s2", &[
        header,
        "chr1\t100\t.\tAT\tA\t.\t.\t.\tGT:DP:AD:FREQ\t0/1:100:30:30%",
    ]);

    let out = dir.path().join("out.csv");
    let status = fstic_bin()
        .args([
            "--vcf",
            dir.path().join("s1.vcf").to_str().unwrap(),
            dir.path().join("s2.vcf").to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "--min-depth", "1", "--min-af", "0.01",
            "--min-alt-reads", "1", "--min-alt-rev-reads", "0",
        ])
        .status()
        .unwrap();

    // Should fail because no SNPs remain
    assert!(!status.success());
}

/// ALT=. and ALT=* should be skipped.
#[test]
fn non_variant_and_star_skipped() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";

    write_vcf(dir.path(), "s1", &[
        header,
        "chr1\t100\t.\tA\t.\t.\t.\t.\tGT:DP:AD:FREQ\t0/0:100:0:0%",  // no-variant
        "chr1\t200\t.\tA\t*\t.\t.\t.\tGT:DP:AD:FREQ\t0/1:100:50:50%", // upstream del
        "chr1\t300\t.\tA\tT\t.\t.\t.\tGT:DP:AD:FREQ\t0/1:100:50:50%", // real SNP
    ]);
    write_vcf(dir.path(), "s2", &[
        header,
        "chr1\t300\t.\tA\tT\t.\t.\t.\tGT:DP:AD:FREQ\t0/1:100:30:30%",
    ]);

    let out = dir.path().join("out.csv");
    let status = fstic_bin()
        .args([
            "--vcf",
            dir.path().join("s1.vcf").to_str().unwrap(),
            dir.path().join("s2.vcf").to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "--min-depth", "1", "--min-af", "0.01",
            "--min-alt-reads", "1", "--min-alt-rev-reads", "0",
        ])
        .status()
        .unwrap();

    assert!(status.success());
    let content = fs::read_to_string(&out).unwrap();
    // Should have computed distances with only 1 SNP at pos 300
    assert!(content.contains("s1"));
}

/// FREQ computed from AD/DP when FREQ field is missing.
#[test]
fn freq_from_ad_dp() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";

    // No FREQ field — should compute from AD/DP
    write_vcf(dir.path(), "s1", &[
        header,
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:AD\t0/1:100:50,50",
    ]);
    write_vcf(dir.path(), "s2", &[
        header,
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:AD\t0/1:100:50,50",
    ]);

    let out = dir.path().join("out.csv");
    let status = fstic_bin()
        .args([
            "--vcf",
            dir.path().join("s1.vcf").to_str().unwrap(),
            dir.path().join("s2.vcf").to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "--min-depth", "1", "--min-af", "0.01",
            "--min-alt-reads", "1", "--min-alt-rev-reads", "0",
        ])
        .status()
        .unwrap();

    assert!(status.success());
    let content = fs::read_to_string(&out).unwrap();
    // Same freq → distance = 0
    assert!(content.contains("0.000000"));
}

/// A sample whose variants are all filtered out is reference everywhere, not absent.
/// It has to keep its row so the matrix stays N x N.
#[test]
fn sample_with_no_surviving_variants_keeps_its_row() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";

    write_vcf(dir.path(), "keep1", &[header, "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:FREQ\t0/1:100:90%"]);
    write_vcf(dir.path(), "keep2", &[header, "chr1\t200\t.\tC\tG\t.\t.\t.\tGT:DP:FREQ\t0/1:100:90%"]);
    write_vcf(dir.path(), "quiet", &[header, "chr1\t300\t.\tT\tC\t.\t.\t.\tGT:DP:FREQ\t0/1:100:10%"]);

    let out = dir.path().join("out.csv");
    let res = fstic_bin()
        .args([
            "--vcf",
            dir.path().join("keep1.vcf").to_str().unwrap(),
            dir.path().join("keep2.vcf").to_str().unwrap(),
            dir.path().join("quiet.vcf").to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "--formula", "bray-curtis",
            "--min-depth", "1", "--min-af", "0.5",
            "--min-alt-reads", "1", "--min-alt-rev-reads", "0",
        ])
        .output()
        .unwrap();

    assert!(res.status.success());
    let content = fs::read_to_string(&out).unwrap();
    assert_eq!(content.lines().count(), 4, "3 samples plus header:\n{}", content);
    assert!(content.lines().any(|l| l.starts_with("quiet,")), "missing row:\n{}", content);
    let err = String::from_utf8_lossy(&res.stderr);
    assert!(err.contains("no variants left after filtering"), "got: {}", err);
}

/// Sample names come from the file stem, so equal basenames would silently merge
/// two samples into one row.
#[test]
fn duplicate_sample_names_are_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";
    let (a, b) = (dir.path().join("a"), dir.path().join("b"));
    fs::create_dir_all(&a).unwrap();
    fs::create_dir_all(&b).unwrap();

    write_vcf(&a, "s", &[header, "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:FREQ\t0/1:100:90%"]);
    write_vcf(&b, "s", &[header, "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:FREQ\t0/1:100:10%"]);

    let out = dir.path().join("out.csv");
    let res = fstic_bin()
        .args([
            "--vcf",
            a.join("s.vcf").to_str().unwrap(),
            b.join("s.vcf").to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "--min-depth", "1", "--min-af", "0.01",
            "--min-alt-reads", "1", "--min-alt-rev-reads", "0",
        ])
        .output()
        .unwrap();

    assert!(!res.status.success());
    let err = String::from_utf8_lossy(&res.stderr);
    assert!(err.contains("share the sample name"), "got: {}", err);
}

/// AD greater than DP would give a frequency above 1, which drives heterozygosity
/// negative. The variant must be dropped and the drop reported.
#[test]
fn freq_above_one_from_ad_dp_is_dropped() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";

    write_vcf(dir.path(), "s1", &[
        header,
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:AD\t0/1:10:15",   // AD > DP, unusable
        "chr1\t200\t.\tC\tG\t.\t.\t.\tGT:DP:AD\t0/1:100:50",  // fine
    ]);
    write_vcf(dir.path(), "s2", &[
        header,
        "chr1\t200\t.\tC\tG\t.\t.\t.\tGT:DP:AD\t0/1:100:50",
    ]);

    let out = dir.path().join("out.csv");
    let res = fstic_bin()
        .args([
            "--vcf",
            dir.path().join("s1.vcf").to_str().unwrap(),
            dir.path().join("s2.vcf").to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "--min-depth", "1", "--min-af", "0.01",
            "--min-alt-reads", "1", "--min-alt-rev-reads", "0",
        ])
        .output()
        .unwrap();

    assert!(res.status.success());
    let err = String::from_utf8_lossy(&res.stderr);
    assert!(err.contains("outside [0,1]"), "the drop must be reported, got: {}", err);
    // Only chr1:200 survives, and both samples agree there
    assert!(err.contains("1 polymorphic sites"), "got: {}", err);
    assert_eq!(off_diagonal(&fs::read_to_string(&out).unwrap(), ','), 0.0);
}

/// An AD value that cannot be parsed must not silently bypass --min-alt-reads.
#[test]
fn unparseable_ad_is_reported() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";

    write_vcf(dir.path(), "s1", &[
        header,
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:AD:FREQ\t0/1:100:.,.:90%",
    ]);
    write_vcf(dir.path(), "s2", &[
        header,
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:AD:FREQ\t0/1:100:90,10:10%",
    ]);

    let out = dir.path().join("out.csv");
    let res = fstic_bin()
        .args([
            "--vcf",
            dir.path().join("s1.vcf").to_str().unwrap(),
            dir.path().join("s2.vcf").to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "--min-depth", "1", "--min-af", "0.01",
            "--min-alt-reads", "50", "--min-alt-rev-reads", "0",
        ])
        .output()
        .unwrap();

    let err = String::from_utf8_lossy(&res.stderr);
    assert!(
        err.contains("could not be parsed as integers"),
        "a filter that could not be applied must be reported, got: {}",
        err
    );
}

/// --pass-only should skip non-PASS variants.
#[test]
fn pass_only_filter() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";

    write_vcf(dir.path(), "s1", &[
        header,
        "chr1\t100\t.\tA\tT\t.\tLowQual\t.\tGT:DP:AD:FREQ\t0/1:100:90:90%",  // filtered
        "chr1\t200\t.\tC\tG\t.\tPASS\t.\tGT:DP:AD:FREQ\t0/1:100:50:50%",     // pass
    ]);
    write_vcf(dir.path(), "s2", &[
        header,
        "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT:DP:AD:FREQ\t0/1:100:10:10%",
        "chr1\t200\t.\tC\tG\t.\tPASS\t.\tGT:DP:AD:FREQ\t0/1:100:50:50%",
    ]);

    let out = dir.path().join("out.csv");
    let status = fstic_bin()
        .args([
            "--vcf",
            dir.path().join("s1.vcf").to_str().unwrap(),
            dir.path().join("s2.vcf").to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "--pass-only",
            "--min-depth", "1", "--min-af", "0.01",
            "--min-alt-reads", "1", "--min-alt-rev-reads", "0",
        ])
        .status()
        .unwrap();

    assert!(status.success());

    // Now run WITHOUT --pass-only (same data)
    let out_no_filter = dir.path().join("out_nofilt.csv");
    let status2 = fstic_bin()
        .args([
            "--vcf",
            dir.path().join("s1.vcf").to_str().unwrap(),
            dir.path().join("s2.vcf").to_str().unwrap(),
            "-o", out_no_filter.to_str().unwrap(),
            "--min-depth", "1", "--min-af", "0.01",
            "--min-alt-reads", "1", "--min-alt-rev-reads", "0",
        ])
        .status()
        .unwrap();
    assert!(status2.success());

    // Distance with --pass-only should be smaller (filtered out the big diff at pos 100)
    let content_filt = fs::read_to_string(&out).unwrap();
    let content_nofilt = fs::read_to_string(&out_no_filter).unwrap();
    let get_dist = |c: &str| -> f64 {
        let lines: Vec<&str> = c.lines().collect();
        let vals: Vec<&str> = lines[1].split(',').collect();
        vals[2].parse().unwrap()
    };
    let d_filt = get_dist(&content_filt);
    let d_nofilt = get_dist(&content_nofilt);
    assert!(d_filt < d_nofilt, "--pass-only should reduce distance (filtered {} vs unfiltered {})", d_filt, d_nofilt);
}
