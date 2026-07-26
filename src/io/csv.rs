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
    let error_count = AtomicUsize::new(0);

    let all_rows: Vec<TableInputRow> = files
        .par_iter()
        .flat_map(|file| {
            let delimiter = match file.extension().and_then(OsStr::to_str) {
                Some("tsv") | Some("tab") => b'\t',
                _ => b',',
            };

            let rdr = csv::ReaderBuilder::new()
                .flexible(false)
                .delimiter(delimiter)
                .from_path(file);

            match rdr {
                Ok(mut reader) => {
                    let mut rows = Vec::new();
                    for result in reader.deserialize() {
                        match result {
                            Ok(row) => rows.push(row),
                            Err(_) => {
                                error_count.fetch_add(1, Ordering::Relaxed);
                            }
                        }
                    }
                    rows
                }
                Err(e) => {
                    eprintln!("Error: cannot open table file {}: {}", file.display(), e);
                    Vec::new()
                }
            }
        })
        .collect();

    let skipped = error_count.load(Ordering::Relaxed);
    if skipped > 0 {
        eprintln!(
            "Warning: skipped {} malformed row(s) in table input.",
            skipped
        );
    }

    let filtered_rows: Vec<TableInputRow> = all_rows.into_iter().filter(|row| {
        // Validate frequency range
        if !row.frequency.is_finite() || !(0.0..=1.0).contains(&row.frequency) {
            return false;
        }
        let depth_ok = row.total_dp.is_none_or(|dp| dp >= criteria.min_depth);
        let freq_ok = row.frequency >= criteria.min_freq;
        let alt_reads_ok = row.alt_dp.is_none_or(|ad| ad >= criteria.min_alt_reads);
        let alt_rev_reads_ok = row
            .alt_rv
            .is_none_or(|arv| arv >= criteria.min_alt_rev_reads);
        depth_ok && freq_ok && alt_reads_ok && alt_rev_reads_ok
    }).collect();

    let mut variants_by_sample: SampleVariants = HashMap::new();
    let mut all_positions: HashSet<GenomicPos> = HashSet::new();

    for row in filtered_rows {
        let gpos = GenomicPos {
            chrom: row.chrom,
            pos: row.position,
        };
        all_positions.insert(gpos.clone());
        let sample_map = variants_by_sample.entry(row.sample).or_default();
        let site_data = sample_map.entry(gpos).or_default();

        if let Some(ref_a) = row.ref_allele {
            site_data.reference_allele = ref_a;
        }
        site_data.freqs.insert(row.sequence, row.frequency);
    }
    Ok((all_positions, variants_by_sample))
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
