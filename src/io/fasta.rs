use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Error, ErrorKind};

/// Reference genome: contig name → sequence bytes (sorted by name for determinism).
pub type ReferenceGenome = BTreeMap<String, Vec<u8>>;

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
                genome.insert(std::mem::take(&mut current_name), std::mem::take(&mut current_seq));
            }
            current_name = header.split_whitespace().next().unwrap_or("").to_string();
        } else {
            // Uppercase to handle soft-masked references
            current_seq.extend(trimmed.bytes().map(|b| b.to_ascii_uppercase()));
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
