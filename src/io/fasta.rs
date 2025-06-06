use std::fs::File;
use std::io::{BufRead, BufReader, Error};

pub fn read_reference_sequence(path: &str) -> Result<Vec<char>, Error> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut sequence = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if !line.starts_with('>') {
            sequence.extend(line.trim().chars());
        }
    }
    Ok(sequence)
}
