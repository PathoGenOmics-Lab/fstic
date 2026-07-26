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
    let mut current_seq: Vec<u8> = Vec::new();
    let mut record_open = false;

    // Coordinates are read straight off these sequences, so a record that silently
    // absorbs someone else's bases shifts every position after it. Reject the
    // malformed shapes instead of guessing.
    let flush = |genome: &mut ReferenceGenome, name: String, seq: Vec<u8>| -> Result<(), Error> {
        if genome.contains_key(&name) {
            return Err(Error::new(
                ErrorKind::InvalidData,
                format!("Duplicate contig name '{}' in FASTA file '{}'", name, path),
            ));
        }
        genome.insert(name, seq);
        Ok(())
    };

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim().trim_start_matches('\u{feff}');
        if trimmed.is_empty() {
            continue;
        }
        if let Some(header) = trimmed.strip_prefix('>') {
            let name = header.split_whitespace().next().unwrap_or("");
            if name.is_empty() {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    format!("FASTA record with no identifier in '{}'", path),
                ));
            }
            if record_open {
                flush(&mut genome, std::mem::take(&mut current_name), std::mem::take(&mut current_seq))?;
            }
            current_seq.clear();
            current_name = name.to_string();
            record_open = true;
        } else {
            if !record_open {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    format!("Sequence data before the first '>' header in FASTA file '{}'", path),
                ));
            }
            // Uppercase to handle soft-masked references
            current_seq.extend(trimmed.bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    if record_open {
        flush(&mut genome, current_name, current_seq)?;
    }

    if genome.is_empty() {
        return Err(Error::new(ErrorKind::InvalidData, "No sequences found in FASTA file"));
    }

    Ok(genome)
}
