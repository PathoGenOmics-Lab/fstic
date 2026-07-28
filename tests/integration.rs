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
    // The diagonal is hard-coded to 0, so a substring check would pass on any
    // matrix at all. Read the off-diagonal cell.
    assert_eq!(off_diagonal(&content, ','), 0.0, "matrix was:\n{}", content);
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
    let res = fstic_bin()
        .args([
            "--vcf",
            dir.path().join("s1.vcf").to_str().unwrap(),
            dir.path().join("s2.vcf").to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "--formula", "bray-curtis",
            "--min-depth", "1", "--min-af", "0.01",
            "--min-alt-reads", "1", "--min-alt-rev-reads", "0",
        ])
        .output()
        .unwrap();

    assert!(res.status.success());
    // Exactly one usable site, chr1:300. Asserting the sample name appears would be
    // satisfied by the header alone.
    let err = String::from_utf8_lossy(&res.stderr);
    assert!(err.contains("1 polymorphic sites"), "got: {}", err);
    // 0.5 against 0.3 over one locus
    let content = fs::read_to_string(&out).unwrap();
    assert!((off_diagonal(&content, ',') - 0.2).abs() < 1e-12, "matrix was:\n{}", content);
}

/// Nei and Reynolds return infinity for a fixed difference, which has to reach the
/// file as NA rather than "inf" for R and pandas to read it.
#[test]
fn non_finite_distance_is_written_as_na() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";

    write_vcf(dir.path(), "fixA", &[header, "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:FREQ\t1/1:100:100%"]);
    write_vcf(dir.path(), "fixB", &[header, "chr1\t100\t.\tA\tG\t.\t.\t.\tGT:DP:FREQ\t1/1:100:100%"]);

    let out = dir.path().join("out.csv");
    let status = fstic_bin()
        .args([
            "--vcf",
            dir.path().join("fixA.vcf").to_str().unwrap(),
            dir.path().join("fixB.vcf").to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "--formula", "nei",
            "--min-depth", "1", "--min-af", "0.01",
            "--min-alt-reads", "1", "--min-alt-rev-reads", "0",
        ])
        .status()
        .unwrap();

    assert!(status.success());
    let content = fs::read_to_string(&out).unwrap();
    assert!(content.contains("NA"), "expected NA, got:\n{}", content);
    assert!(!content.to_lowercase().contains("inf"), "got:\n{}", content);
}

/// A sample name containing the output separator has to be quoted, or the matrix
/// gains a column and every downstream parser misreads it.
#[test]
fn sample_names_containing_the_separator_are_quoted() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";

    write_vcf(dir.path(), "with,comma", &[header, "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:FREQ\t0/1:100:90%"]);
    write_vcf(dir.path(), "plain", &[header, "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:FREQ\t0/1:100:10%"]);

    let out = dir.path().join("out.csv");
    let status = fstic_bin()
        .args([
            "--vcf",
            dir.path().join("with,comma.vcf").to_str().unwrap(),
            dir.path().join("plain.vcf").to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "--min-depth", "1", "--min-af", "0.01",
            "--min-alt-reads", "1", "--min-alt-rev-reads", "0",
        ])
        .status()
        .unwrap();

    assert!(status.success());
    let content = fs::read_to_string(&out).unwrap();
    assert!(content.contains("\"with,comma\""), "name was not quoted:\n{}", content);
    // Header plus two rows, each with 3 fields once quoting is honoured
    for line in content.lines() {
        let fields = line.matches(',').count() - line.matches("\"with,comma\"").count();
        assert_eq!(fields, 2, "row has the wrong field count: {}", line);
    }
}

/// FREQ computed from AD/DP when the FREQ field is missing.
///
/// AD is "ref,alt", so the second component is the one to use. A symmetric fixture
/// like AD=50,50 in both samples cannot catch a ref/alt swap, and neither can two
/// samples with mirrored values, since the distance is the same either way. Pairing
/// an AD-derived frequency against an explicit FREQ pins the absolute value.
#[test]
fn freq_from_ad_dp() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";

    // AD=20,80 means 80 alt reads out of 100, i.e. 0.8
    write_vcf(dir.path(), "s1", &[
        header,
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:AD\t0/1:100:20,80",
    ]);
    // Same frequency, stated outright
    write_vcf(dir.path(), "s2", &[
        header,
        "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:FREQ\t0/1:100:80%",
    ]);

    let out = dir.path().join("out.csv");
    let status = fstic_bin()
        .args([
            "--vcf",
            dir.path().join("s1.vcf").to_str().unwrap(),
            dir.path().join("s2.vcf").to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "--formula", "bray-curtis",
            "--min-depth", "1", "--min-af", "0.01",
            "--min-alt-reads", "1", "--min-alt-rev-reads", "0",
        ])
        .status()
        .unwrap();

    assert!(status.success());
    let content = fs::read_to_string(&out).unwrap();
    // Reading AD as the ref depth would give 0.2 against 0.8, i.e. 0.6 here.
    assert_eq!(off_diagonal(&content, ','), 0.0, "matrix was:\n{}", content);
}

fn write_ref(dir: &Path) -> std::path::PathBuf {
    let p = dir.join("ref.fa");
    fs::write(&p, format!(">chr1\n{}\n", "A".repeat(500))).unwrap();
    p
}

fn run_table(dir: &Path, body: &str, extra: &[&str]) -> (std::process::Output, std::path::PathBuf) {
    let table = dir.join("t.csv");
    fs::write(&table, body).unwrap();
    let reference = write_ref(dir);
    let out = dir.join("out.csv");
    let mut args: Vec<String> = vec![
        "--table".into(), table.to_str().unwrap().into(),
        "-r".into(), reference.to_str().unwrap().into(),
        "-o".into(), out.to_str().unwrap().into(),
        "--min-depth".into(), "1".into(),
        "--min-af".into(), "0.01".into(),
        "--min-alt-reads".into(), "1".into(),
        "--min-alt-rev-reads".into(), "0".into(),
    ];
    args.extend(extra.iter().map(|s| s.to_string()));
    (fstic_bin().args(&args).output().unwrap(), out)
}

/// Table mode must agree with VCF mode on the same underlying data.
#[test]
fn table_matches_vcf_on_equivalent_data() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";
    write_vcf(dir.path(), "sA", &[header, "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:FREQ\t0/1:100:50%"]);
    write_vcf(dir.path(), "sB", &[header, "chr1\t100\t.\tA\tT\t.\t.\t.\tGT:DP:FREQ\t0/1:100:80%"]);

    let vcf_out = dir.path().join("vcf.csv");
    let status = fstic_bin()
        .args([
            "--vcf",
            dir.path().join("sA.vcf").to_str().unwrap(),
            dir.path().join("sB.vcf").to_str().unwrap(),
            "-o", vcf_out.to_str().unwrap(),
            "--formula", "chord",
            "--min-depth", "1", "--min-af", "0.01",
            "--min-alt-reads", "1", "--min-alt-rev-reads", "0",
        ])
        .status()
        .unwrap();
    assert!(status.success());

    let (res, table_out) = run_table(
        dir.path(),
        "sample,chrom,position,ref_allele,sequence,frequency\n\
         sA,chr1,100,A,T,0.5\n\
         sB,chr1,100,A,T,0.8\n",
        &["--formula", "chord"],
    );
    assert!(res.status.success(), "{}", String::from_utf8_lossy(&res.stderr));

    let from_vcf = off_diagonal(&fs::read_to_string(&vcf_out).unwrap(), ',');
    let from_table = off_diagonal(&fs::read_to_string(&table_out).unwrap(), ',');
    assert!((from_vcf - from_table).abs() < 1e-12, "{} vs {}", from_vcf, from_table);
}

/// The README promises case-insensitive column names.
#[test]
fn table_headers_are_case_insensitive() {
    let dir = tempfile::tempdir().unwrap();
    let (res, out) = run_table(
        dir.path(),
        "Sample,Position,Ref_Allele,Sequence,Frequency\n\
         sA,100,A,T,0.5\n\
         sB,100,A,T,0.8\n",
        &["--formula", "bray-curtis"],
    );
    assert!(res.status.success(), "{}", String::from_utf8_lossy(&res.stderr));
    assert!((off_diagonal(&fs::read_to_string(&out).unwrap(), ',') - 0.3).abs() < 1e-12);
}

/// Percentages are auto-detected in the VCF FREQ field, so tables must match.
#[test]
fn table_accepts_percentage_frequencies() {
    let dir = tempfile::tempdir().unwrap();
    let (res, out) = run_table(
        dir.path(),
        "sample,position,ref_allele,sequence,frequency\n\
         sA,100,A,T,50%\n\
         sB,100,A,T,80%\n",
        &["--formula", "bray-curtis"],
    );
    assert!(res.status.success(), "{}", String::from_utf8_lossy(&res.stderr));
    assert!((off_diagonal(&fs::read_to_string(&out).unwrap(), ',') - 0.3).abs() < 1e-12);
}

/// Contradictory duplicate rows keep the first value, like the VCF reader, and say so.
#[test]
fn table_duplicate_rows_keep_first_and_warn() {
    let dir = tempfile::tempdir().unwrap();
    let (res, out) = run_table(
        dir.path(),
        "sample,position,ref_allele,sequence,frequency\n\
         sA,100,A,T,0.5\n\
         sA,100,A,T,0.9\n\
         sB,100,A,T,0.8\n",
        &["--formula", "bray-curtis"],
    );
    assert!(res.status.success());
    let err = String::from_utf8_lossy(&res.stderr);
    assert!(err.contains("duplicate"), "got: {}", err);
    // 0.5 kept, not 0.9
    assert!((off_diagonal(&fs::read_to_string(&out).unwrap(), ',') - 0.3).abs() < 1e-12);
}

/// Mixing tables with and without a chrom column would split one site into two loci.
#[test]
fn table_mixed_chrom_column_is_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let reference = write_ref(dir.path());
    let a = dir.path().join("a.csv");
    let b = dir.path().join("b.csv");
    fs::write(&a, "sample,chrom,position,ref_allele,sequence,frequency\nsA,chr1,100,A,T,0.5\n").unwrap();
    fs::write(&b, "sample,position,ref_allele,sequence,frequency\nsB,100,A,T,0.8\n").unwrap();

    let out = dir.path().join("out.csv");
    let res = fstic_bin()
        .args([
            "--table", a.to_str().unwrap(), b.to_str().unwrap(),
            "-r", reference.to_str().unwrap(),
            "-o", out.to_str().unwrap(),
            "--min-depth", "1", "--min-af", "0.01",
            "--min-alt-reads", "1", "--min-alt-rev-reads", "0",
        ])
        .output()
        .unwrap();

    assert!(!res.status.success());
    let err = String::from_utf8_lossy(&res.stderr);
    assert!(err.contains("chrom"), "got: {}", err);
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

/// Output must not depend on the thread count: a 4-core and a 16-core machine
/// running the same command have to produce the same file.
#[test]
fn output_is_independent_of_worker_count() {
    let dir = tempfile::tempdir().unwrap();
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample";

    // Enough loci, at irregular frequencies, for the summation order to show up.
    let mut a = vec![header.to_string()];
    let mut b = vec![header.to_string()];
    for pos in 1..4000 {
        let fa = 3.0 + (pos as f64 * 7.3) % 91.0;
        let fb = 2.0 + (pos as f64 * 11.7) % 93.0;
        a.push(format!("chr1\t{}\t.\tA\tT\t.\t.\t.\tGT:DP:FREQ\t0/1:100:{:.6}%", pos, fa));
        b.push(format!("chr1\t{}\t.\tA\tT\t.\t.\t.\tGT:DP:FREQ\t0/1:100:{:.6}%", pos, fb));
    }
    let a_ref: Vec<&str> = a.iter().map(|s| s.as_str()).collect();
    let b_ref: Vec<&str> = b.iter().map(|s| s.as_str()).collect();
    write_vcf(dir.path(), "wa", &a_ref);
    write_vcf(dir.path(), "wb", &b_ref);

    let run = |workers: &str, name: &str| -> String {
        let out = dir.path().join(name);
        let status = fstic_bin()
            .args([
                "--vcf",
                dir.path().join("wa.vcf").to_str().unwrap(),
                dir.path().join("wb.vcf").to_str().unwrap(),
                "-o", out.to_str().unwrap(),
                "--formula", "fst",
                "--workers", workers,
                "--min-depth", "1", "--min-af", "0.01",
                "--min-alt-reads", "1", "--min-alt-rev-reads", "0",
            ])
            .status()
            .unwrap();
        assert!(status.success());
        fs::read_to_string(&out).unwrap()
    };

    let one = run("1", "w1.csv");
    for (workers, name) in [("2", "w2.csv"), ("3", "w3.csv"), ("8", "w8.csv")] {
        assert_eq!(one, run(workers, name), "--workers {} changed the output", workers);
    }
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
