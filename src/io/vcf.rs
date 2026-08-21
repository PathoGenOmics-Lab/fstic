use crate::types::{FilterCriteria, GenomicPos, SampleVariants, VcfVariant};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

/// Counters for input that was dropped or could not be used, reported once per run
/// rather than once per line: a whole-genome VCF with a systematic problem would
/// otherwise bury the terminal in millions of identical warnings.
#[derive(Default)]
pub struct VcfStats {
    multi_allelic: AtomicUsize,
    bad_int: AtomicUsize,
    bad_freq: AtomicUsize,
    freq_out_of_range: AtomicUsize,
}

impl VcfStats {
    fn report(&self) {
        let n = self.multi_allelic.load(Ordering::Relaxed);
        if n > 0 {
            eprintln!(
                "Warning: skipped {} multi-allelic site(s) (>1 ALT allele). \
                 Only bi-allelic sites are currently supported.",
                n
            );
        }
        let n = self.bad_int.load(Ordering::Relaxed);
        if n > 0 {
            eprintln!(
                "Warning: {} FORMAT value(s) for DP/AD/ADR could not be parsed as integers. \
                 The corresponding depth filters could not be applied to those variants.",
                n
            );
        }
        let n = self.bad_freq.load(Ordering::Relaxed);
        if n > 0 {
            eprintln!(
                "Warning: {} variant(s) had an unparseable FREQ value and no usable AD/DP \
                 fallback; they were dropped.",
                n
            );
        }
        let n = self.freq_out_of_range.load(Ordering::Relaxed);
        if n > 0 {
            eprintln!(
                "Warning: dropped {} variant(s) whose allele frequency fell outside [0,1] \
                 (a FREQ out of range, or AD greater than DP).",
                n
            );
        }
    }
}

pub fn read_vcf_files(
    files: &[PathBuf],
    criteria: &FilterCriteria,
) -> Result<(HashSet<GenomicPos>, SampleVariants), String> {
    let stats = VcfStats::default();

    let parsed: Vec<ParsedVcf> = files
        .par_iter()
        .map(|file| parse_and_filter_vcf(file, criteria, &stats))
        .collect();

    stats.report();

    // The sample name is the file stem, so two inputs with the same basename in
    // different directories would land in the same matrix row and one would be lost.
    let mut seen: HashMap<&str, &PathBuf> = HashMap::new();
    for p in &parsed {
        if let Some(first) = seen.insert(p.sample.as_str(), &p.path) {
            return Err(format!(
                "Two input files share the sample name '{}': {} and {}. \
                 Sample names come from the file name, so rename one of them.",
                p.sample,
                first.display(),
                p.path.display()
            ));
        }
    }

    for p in &parsed {
        if let Some(err) = &p.io_error {
            return Err(format!(
                "Error reading VCF file {}: {}",
                p.path.display(),
                err
            ));
        }
        if p.data_lines == 0 {
            return Err(format!(
                "No VCF data lines could be read from {}. \
                 Check the file is an uncompressed, tab-delimited VCF.",
                p.path.display()
            ));
        }
        if p.short_lines > 0 {
            eprintln!(
                "Warning: {}: {} data line(s) had fewer than 10 tab-separated columns \
                 and were skipped (sites-only or space-delimited VCF?).",
                p.path.display(),
                p.short_lines
            );
        }
        if p.variants.is_empty() {
            eprintln!(
                "Warning: {} has no variants left after filtering; \
                 '{}' is treated as identical to the reference.",
                p.path.display(),
                p.sample
            );
        }
    }

    let mut variants_by_sample: SampleVariants = HashMap::new();
    let mut all_positions: HashSet<GenomicPos> = HashSet::new();

    // Register every input up front so a sample with no surviving variants still
    // gets a matrix row rather than vanishing from the output.
    for p in &parsed {
        variants_by_sample.entry(p.sample.clone()).or_default();
    }

    for variant in parsed.into_iter().flat_map(|p| p.variants) {
        let gpos = GenomicPos {
            chrom: variant.chrom.clone(),
            pos: variant.pos,
        };
        all_positions.insert(gpos.clone());
        let sample_map = variants_by_sample
            .entry(variant.sample.clone())
            .or_default();

        let site_data = sample_map.entry(gpos).or_default();
        // Only set ref_allele if not already set (avoid overwriting from duplicate lines)
        if site_data.reference_allele.is_empty() {
            site_data.reference_allele = variant.ref_allele;
        }
        // For duplicate positions: first alt wins (second line at same pos is likely
        // an artifact from split multi-allelics or overlapping indel/SNP calls)
        site_data
            .freqs
            .entry(variant.alt_allele)
            .or_insert(variant.alt_freq);
    }

    Ok((all_positions, variants_by_sample))
}

/// One input file's worth of variants, plus what happened while reading it.
struct ParsedVcf {
    sample: String,
    path: PathBuf,
    variants: Vec<VcfVariant>,
    /// Non-header lines seen. Zero means the file was not a readable plain VCF.
    data_lines: usize,
    short_lines: usize,
    io_error: Option<String>,
}

fn parse_and_filter_vcf(file: &Path, criteria: &FilterCriteria, stats: &VcfStats) -> ParsedVcf {
    // Sample ID from filename (one VCF per sample is the expected input model)
    let sample_id = file
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown_sample")
        .to_string();

    let mut report = ParsedVcf {
        sample: sample_id.clone(),
        path: file.to_path_buf(),
        variants: Vec::new(),
        data_lines: 0,
        short_lines: 0,
        io_error: None,
    };

    let vcf_file = match File::open(file) {
        Ok(f) => f,
        Err(e) => {
            report.io_error = Some(e.to_string());
            return report;
        }
    };
    let reader = BufReader::new(vcf_file);
    let mut warned: HashSet<String> = HashSet::new();

    let mut variants = Vec::new();
    for line in reader.lines() {
        // Bailing out quietly here truncated the file at the first read error, and a
        // gzipped VCF hits one immediately: the sample came out empty, exit code 0.
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                report.io_error = Some(e.to_string());
                return report;
            }
        };
        if line.starts_with("##") {
            continue;
        }

        // Parse #CHROM header to check for multi-sample VCFs
        if line.starts_with("#CHROM") {
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() > 10 {
                eprintln!(
                    "Warning: VCF file '{}' contains {} samples. \
                     Only the first sample column will be used. \
                     Split into per-sample VCFs for full analysis.",
                    file.display(),
                    cols.len() - 9
                );
            }
            continue;
        }

        report.data_lines += 1;

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            report.short_lines += 1;
            continue;
        }

        // FILTER column check (index 6)
        if criteria.pass_only {
            let filter = fields[6];
            if filter != "PASS" && filter != "." {
                continue;
            }
        }

        // Skip non-variant sites (ALT = ".")
        if fields[4] == "." {
            continue;
        }

        // Skip upstream deletion markers (ALT = "*")
        if fields[4] == "*" || fields[4].contains('*') {
            continue;
        }

        // Skip multi-allelic sites with a count
        if fields[4].contains(',') {
            stats.multi_allelic.fetch_add(1, Ordering::Relaxed);
            continue;
        }

        // Skip indels: only process SNPs (single-base REF and ALT)
        if fields[3].len() != 1 || fields[4].len() != 1 {
            continue;
        }

        if let Some(variant) = parse_vcf_line(&fields, &sample_id, file, &mut warned, stats) {
            variants.push(variant);
        }
    }

    report.variants = variants
        .into_iter()
        .filter(|v| {
            let depth_ok = v.total_dp.is_none_or(|dp| dp >= criteria.min_depth);
            let freq_ok = v.alt_freq >= criteria.min_freq;
            let alt_reads_ok = v.alt_dp.is_none_or(|ad| ad >= criteria.min_alt_reads);
            let alt_rev_reads_ok = v.alt_rv.is_none_or(|arv| arv >= criteria.min_alt_rev_reads);
            depth_ok && freq_ok && alt_reads_ok && alt_rev_reads_ok
        })
        .collect();
    report
}

fn parse_vcf_line(
    fields: &[&str],
    sample_id: &str,
    filepath: &Path,
    warned: &mut HashSet<String>,
    stats: &VcfStats,
) -> Option<VcfVariant> {
    let chrom = fields[0].to_string();
    let pos = fields[1].parse::<usize>().ok()?;
    let ref_allele = fields[3].to_string();

    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let sample_values: Vec<&str> = fields[9].split(':').collect();
    let format_map: HashMap<&str, &str> = format_keys
        .iter()
        .copied()
        .zip(sample_values.iter().copied())
        .collect();

    let get_value = |key: &str| -> &str { format_map.get(key).copied().unwrap_or(".") };

    // A value that is present but unparseable is not the same as an absent field: the
    // filter the user asked for cannot be applied either way, but corrupt input has to
    // be counted so the run does not look clean.
    let parse_int_warn = |key: &str, warned: &mut HashSet<String>| -> Option<u32> {
        match get_value(key) {
            "." => {
                if warned.insert(key.to_string()) {
                    eprintln!(
                        "Warning: VCF file '{}' is missing FORMAT field '{}'. \
                         Filtering on this field will be skipped.",
                        filepath.display(),
                        key
                    );
                }
                None
            }
            val => {
                let parsed = val.parse::<u32>().ok();
                if parsed.is_none() {
                    stats.bad_int.fetch_add(1, Ordering::Relaxed);
                }
                parsed
            }
        }
    };

    let total_dp = parse_int_warn("DP", warned);

    // AD field: VCF spec uses comma-separated "ref_depth,alt_depth[,alt2_depth,...]"
    // VarScan uses a single integer for alt reads only.
    // Handle both formats.
    let alt_dp = match get_value("AD") {
        "." => {
            if warned.insert("AD".to_string()) {
                eprintln!(
                    "Warning: VCF file '{}' is missing FORMAT field 'AD'. \
                     Filtering on this field will be skipped.",
                    filepath.display(),
                );
            }
            None
        }
        val => {
            let parsed = if val.contains(',') {
                // Standard VCF: AD=ref,alt, take the second value (alt depth)
                val.split(',').nth(1).and_then(|s| s.parse::<u32>().ok())
            } else {
                // VarScan style: AD=alt_reads (single integer)
                val.parse::<u32>().ok()
            };
            if parsed.is_none() {
                stats.bad_int.fetch_add(1, Ordering::Relaxed);
            }
            parsed
        }
    };

    let alt_rv = parse_int_warn("ADR", warned);

    // AD/DP is the fallback whenever FREQ is unusable. AD > DP happens with real
    // callers and would otherwise yield a frequency above 1, which drives
    // heterozygosity negative and takes every metric out of its range.
    let freq_from_depth = || -> Option<f64> {
        match (alt_dp, total_dp) {
            (Some(ad), Some(dp)) if dp > 0 => Some(ad as f64 / dp as f64),
            _ => None,
        }
    };

    // FREQ: handle proportion (0.5), percentage (50%), or compute from AD/DP
    let alt_freq = match get_value("FREQ") {
        "." => {
            if warned.insert("FREQ".to_string()) {
                eprintln!(
                    "Warning: VCF file '{}' is missing FORMAT field 'FREQ'. \
                     Will compute from AD/DP when available.",
                    filepath.display()
                );
            }
            freq_from_depth()
        }
        val => {
            // Number=A fields arrive comma-separated; this is a bi-allelic site by
            // the time we get here, so the first component is the one we want.
            let val = val.split(',').next().unwrap_or(val).trim();
            let parsed = if let Some(stripped) = val.strip_suffix('%') {
                stripped.trim().parse::<f64>().ok().map(|f| f / 100.0)
            } else {
                val.parse::<f64>().ok()
            };
            match parsed {
                Some(f) => Some(f),
                None => {
                    let fallback = freq_from_depth();
                    if fallback.is_none() {
                        stats.bad_freq.fetch_add(1, Ordering::Relaxed);
                    }
                    fallback
                }
            }
        }
    };

    // Whatever the source, an allele frequency outside [0,1] is not usable.
    let alt_freq = match alt_freq {
        Some(f) if f.is_finite() && (0.0..=1.0).contains(&f) => f,
        Some(_) => {
            stats.freq_out_of_range.fetch_add(1, Ordering::Relaxed);
            return None;
        }
        // Cannot compute distances without knowing the allele frequency.
        None => return None,
    };

    Some(VcfVariant {
        sample: sample_id.to_string(),
        chrom,
        pos,
        ref_allele,
        alt_allele: fields[4].to_string(),
        total_dp,
        alt_dp,
        alt_freq, // guaranteed non-None by the ? above
        alt_rv,
    })
}
