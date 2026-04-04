use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Error, ErrorKind};

/// Reference genome: contig name → sequence bytes.
pub type ReferenceGenome = HashMap<String, Vec<u8>>;

/// Reads a multi-contig reference FASTA file into a name→sequence map.
///
/// Stores sequences as `Vec<u8>` (1 byte/base) instead of `Vec<char>` (4 bytes/base).
pub fn read_reference(path: &str) -> Result<ReferenceGenome, Error> {
    let file = File::open(path).map_err(|e| {
        Error::new(ErrorKind::NotFound, format!("Cannot open FASTA file '{}': {}", path, e))
    })?;
    let reader = BufReader::new(file);
    let mut genome = ReferenceGenome::new();
    let mut current_name = String::new();
    let mut current_seq = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(header) = trimmed.strip_prefix('>') {
            if !current_name.is_empty() {
                genome.insert(current_name.clone(), current_seq.clone());
                current_seq.clear();
            }
            current_name = header.split_whitespace().next().unwrap_or("").to_string();
        } else {
            current_seq.extend_from_slice(trimmed.as_bytes());
        }
    }
    if !current_name.is_empty() {
        genome.insert(current_name, current_seq);
    }

    if genome.is_empty() {
        return Err(Error::new(ErrorKind::InvalidData, "No sequences found in FASTA file"));
    }

    Ok(genome)
}

/// Gets a base at 1-based position from the first (or only) contig.
///
/// Convenience for single-contig references used in VCF mode.
pub fn get_base_at(genome: &ReferenceGenome, pos: usize) -> Option<u8> {
    // Try the first contig (most common case: single-contig reference)
    genome.values().next().and_then(|seq| seq.get(pos.saturating_sub(1)).copied())
}
