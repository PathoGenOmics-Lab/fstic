use crate::types::{FilterCriteria, GenomicPos, SampleVariants, TableInputRow};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::ffi::OsStr;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

pub fn read_csv_files(
    files: &[PathBuf],
    criteria: &FilterCriteria,
) -> Result<(HashSet<GenomicPos>, SampleVariants), String> {
    let out_of_range = AtomicUsize::new(0);

    let per_file: Vec<TableFile> = files.par_iter().map(read_one_table).collect();

    let mut all_rows: Vec<TableInputRow> = Vec::new();
    for f in per_file {
        if let Some(err) = f.open_error {
            return Err(format!(
                "Cannot open table file {}: {}",
                f.path.display(),
                err
            ));
        }
        if f.rows.is_empty() && f.error_count > 0 {
            return Err(format!(
                "Every data row in {} failed to parse ({} row(s)). First error: {}",
                f.path.display(),
                f.error_count,
                f.first_error.unwrap_or_else(|| "unknown".to_string())
            ));
        }
        if f.error_count > 0 {
            eprintln!(
                "Warning: skipped {} malformed row(s) in {}. First error: {}",
                f.error_count,
                f.path.display(),
                f.first_error.unwrap_or_else(|| "unknown".to_string())
            );
        }
        all_rows.extend(f.rows);
    }

    // A table with no chrom column and one with a real chrom column cannot be
    // combined: the same physical site would land under two different keys.
    let with_chrom = all_rows.iter().filter(|r| r.chrom.is_some()).count();
    if with_chrom != 0 && with_chrom != all_rows.len() {
        return Err("Some table rows carry a 'chrom' column and others do not. \
             Add the column everywhere, or remove it everywhere."
            .to_string());
    }

    let filtered_rows: Vec<TableInputRow> = all_rows
        .into_iter()
        .filter(|row| {
            // Validate frequency range
            if !row.frequency.is_finite() || !(0.0..=1.0).contains(&row.frequency) {
                out_of_range.fetch_add(1, Ordering::Relaxed);
                return false;
            }
            let depth_ok = row.total_dp.is_none_or(|dp| dp >= criteria.min_depth);
            let freq_ok = row.frequency >= criteria.min_freq;
            let alt_reads_ok = row.alt_dp.is_none_or(|ad| ad >= criteria.min_alt_reads);
            let alt_rev_reads_ok = row
                .alt_rv
                .is_none_or(|arv| arv >= criteria.min_alt_rev_reads);
            depth_ok && freq_ok && alt_reads_ok && alt_rev_reads_ok
        })
        .collect();

    let n = out_of_range.load(Ordering::Relaxed);
    if n > 0 {
        eprintln!(
            "Warning: dropped {} table row(s) whose frequency was outside [0,1] or not finite.",
            n
        );
    }

    let mut variants_by_sample: SampleVariants = HashMap::new();
    let mut all_positions: HashSet<GenomicPos> = HashSet::new();
    let mut duplicates = 0usize;

    for row in filtered_rows {
        let gpos = GenomicPos {
            // Rows either all carry a chrom or none do, so this placeholder is
            // applied uniformly and cannot split a site in two.
            chrom: row
                .chrom
                .filter(|c| !c.is_empty())
                .unwrap_or_else(|| ".".to_string()),
            pos: row.position,
        };
        all_positions.insert(gpos.clone());
        let sample_map = variants_by_sample.entry(row.sample).or_default();
        let site_data = sample_map.entry(gpos).or_default();

        if let Some(ref_a) = row.ref_allele {
            if site_data.reference_allele.is_empty() {
                site_data.reference_allele = ref_a;
            }
        }
        // First value wins, matching the VCF reader. Contradictory duplicates have no
        // defensible answer, so the point is to say they were there.
        if site_data.freqs.contains_key(&row.sequence) {
            duplicates += 1;
        }
        site_data.freqs.entry(row.sequence).or_insert(row.frequency);
    }

    if duplicates > 0 {
        eprintln!(
            "Warning: {} duplicate (sample, chrom, position, allele) row(s); first value kept.",
            duplicates
        );
    }

    Ok((all_positions, variants_by_sample))
}

/// One table file's rows plus what went wrong reading it.
struct TableFile {
    path: PathBuf,
    rows: Vec<TableInputRow>,
    error_count: usize,
    first_error: Option<String>,
    open_error: Option<String>,
}

fn read_one_table(file: &PathBuf) -> TableFile {
    let mut out = TableFile {
        path: file.clone(),
        rows: Vec::new(),
        error_count: 0,
        first_error: None,
        open_error: None,
    };

    let delimiter = match file.extension().and_then(OsStr::to_str) {
        Some("tsv") | Some("tab") => b'\t',
        _ => b',',
    };

    let mut reader = match csv::ReaderBuilder::new()
        .flexible(false)
        .delimiter(delimiter)
        .from_path(file)
    {
        Ok(r) => r,
        Err(e) => {
            out.open_error = Some(e.to_string());
            return out;
        }
    };

    // The README promises case-insensitive column names, and without this a header
    // of "Sample,Position,..." makes every row fail to deserialise, which silently
    // deletes the whole file's worth of samples.
    match reader.headers() {
        Ok(headers) => {
            let lowered: csv::StringRecord = headers
                .iter()
                .map(|h| h.trim().to_ascii_lowercase())
                .collect();
            reader.set_headers(lowered);
        }
        Err(e) => {
            out.open_error = Some(e.to_string());
            return out;
        }
    }

    for result in reader.deserialize() {
        match result {
            Ok(row) => out.rows.push(row),
            Err(e) => {
                out.error_count += 1;
                if out.first_error.is_none() {
                    out.first_error = Some(e.to_string());
                }
            }
        }
    }
    out
}

/// Quote a field if it contains the separator, a quote, or a newline.
fn quote_field(field: &str, sep: &str) -> String {
    if field.contains(sep) || field.contains('"') || field.contains('\n') {
        format!("\"{}\"", field.replace('"', "\"\""))
    } else {
        field.to_string()
    }
}

pub fn write_distance_matrix(
    path: &str,
    matrix: &[Vec<f64>],
    samples: &[String],
    use_tab: bool,
) -> std::io::Result<()> {
    let sep = if use_tab { "\t" } else { "," };
    let mut file = File::create(path)?;

    // Header
    write!(file, "sample")?;
    for s in samples {
        write!(file, "{}{}", sep, quote_field(s, sep))?;
    }
    writeln!(file)?;

    // Rows
    for (i, row) in matrix.iter().enumerate() {
        write!(file, "{}", quote_field(&samples[i], sep))?;
        for val in row {
            if val.is_finite() {
                // 10 decimal places prevents rounding closely-related pairs to 0
                // (e.g. MTB FST ~1e-8). Standard tools accept this fine.
                write!(file, "{}{:.10}", sep, val)?;
            } else {
                // Nei/Reynolds can produce Infinity for fixed differences.
                // Use NA for compatibility with R/Python/downstream tools.
                write!(file, "{}NA", sep)?;
            }
        }
        writeln!(file)?;
    }
    Ok(())
}
