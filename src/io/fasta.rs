use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

/// Multi-contig reference: maps contig name -> sequence bytes.
pub type ReferenceGenome = HashMap<String, Vec<u8>>;

/// Read a reference FASTA file (plain or gzipped) into a multi-contig map.
pub fn read_reference_genome(path: &str) -> Result<ReferenceGenome> {
    let file = File::open(path).with_context(|| format!("Cannot open FASTA file: {}", path))?;
    let reader: Box<dyn Read> = if path.ends_with(".gz") {
        Box::new(flate2::read::GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let buf = BufReader::new(reader);

    let mut genome = ReferenceGenome::new();
    let mut current_contig = String::new();
    let mut current_seq = Vec::new();

    for line in buf.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.starts_with('>') {
            if !current_contig.is_empty() {
                genome.insert(current_contig.clone(), current_seq.clone());
                current_seq.clear();
            }
            current_contig = trimmed[1..].split_whitespace().next().unwrap_or("").to_string();
        } else {
            current_seq.extend(trimmed.as_bytes());
        }
    }
    if !current_contig.is_empty() {
        genome.insert(current_contig, current_seq);
    }

    Ok(genome)
}

/// Legacy single-sequence reader for backwards compatibility.
/// Concatenates all contigs into a single sequence.
pub fn read_reference_sequence(path: &str) -> Result<Vec<char>> {
    let genome = read_reference_genome(path)?;
    let mut sequence = Vec::new();
    // Sort by contig name for deterministic order
    let mut contigs: Vec<_> = genome.into_iter().collect();
    contigs.sort_by(|a, b| a.0.cmp(&b.0));
    for (_, seq) in contigs {
        for &byte in &seq {
            sequence.push(byte as char);
        }
    }
    Ok(sequence)
}
