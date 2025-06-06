use crate::types::{FilterCriteria, SampleVariants, TableInputRow};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::ffi::OsStr;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

pub fn read_csv_files(
    files: &[PathBuf],
    criteria: &FilterCriteria,
) -> (HashSet<usize>, SampleVariants) {
    let all_rows: Vec<TableInputRow> = files
        .par_iter()
        .flat_map(|file| {
            let delimiter = match file.extension().and_then(OsStr::to_str) {
                Some("tsv") | Some("tab") => b'\t',
                _ => b',',
            };

            let mut rdr = csv::ReaderBuilder::new()
                .flexible(true)
                .delimiter(delimiter)
                .from_path(file)
                .expect("Cannot open input table file");

            rdr.deserialize().filter_map(Result::ok).collect::<Vec<_>>()
        })
        .collect();

    // Apply filters to the rows from the table files
    let filtered_rows = all_rows.into_iter().filter(|row| {
        let depth_ok = row.total_dp.map_or(true, |dp| dp >= criteria.min_depth);
        let freq_ok = row.frequency >= criteria.min_freq;
        let alt_reads_ok = row.alt_dp.map_or(true, |ad| ad >= criteria.min_alt_reads);
        let alt_rev_reads_ok = row.alt_rv.map_or(true, |arv| arv >= criteria.min_alt_rev_reads);
        depth_ok && freq_ok && alt_reads_ok && alt_rev_reads_ok
    });

    let mut variants_by_sample: SampleVariants = HashMap::new();
    let mut all_positions: HashSet<usize> = HashSet::new();

    for row in filtered_rows {
        all_positions.insert(row.position);
        let sample_map = variants_by_sample
            .entry(row.sample)
            .or_default();

        let site_data = sample_map
            .entry(row.position)
            .or_default();
        
        if let Some(ref_a) = row.ref_allele {
            site_data.reference_allele = ref_a;
        }
        site_data.freqs.insert(row.sequence, row.frequency);
    }
    (all_positions, variants_by_sample)
}

pub fn write_distance_matrix(
    path: &str,
    matrix: &[Vec<f64>],
    samples: &[String],
) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    writeln!(file, "\t{}", samples.join("\t"))?;

    for (i, row) in matrix.iter().enumerate() {
        let row_str: Vec<String> = row.iter().map(|v| format!("{:.6}", v)).collect();
        writeln!(file, "{}\t{}", samples[i], row_str.join("\t"))?;
    }
    Ok(())
}
