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
