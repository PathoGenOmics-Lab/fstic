use std::fs::File;
use std::io::{BufRead, BufReader, Error};
use std::path::PathBuf;

pub mod vcf;
pub mod csv;
pub mod fasta;

pub fn read_file_list(path: &str) -> Result<Vec<PathBuf>, Error> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut paths = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if !trimmed.is_empty() && !trimmed.starts_with('#') {
            paths.push(PathBuf::from(trimmed));
        }
    }
    Ok(paths)
}
