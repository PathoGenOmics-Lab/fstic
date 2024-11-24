use clap::{App, Arg};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;

/// Struct to hold variant information
#[derive(Clone)]
struct Variant {
    region: String,
    pos: usize,
    ref_allele: String,
    alt_allele: String,
    total_dp: usize,
    ref_dp: usize,
    alt_dp: usize,
    alt_freq: f64,
    alt_rv: usize,
}

fn main() {
    // Parse command-line arguments
    let matches = App::new("Fstigo")
        .version("1.0.0")
        .about("Calculates pairwise Fst values between samples using alelle frequencies")
        .author("Paula Ruiz-Rodriguez <paula.ruiz.rodriguez@csic.es>")
        .arg(
            Arg::with_name("vcf")
                .short("v")
                .long("vcf")
                .value_name("VCF_FILES")
                .help("Input VCF files")
                .multiple(true)
                .conflicts_with("vcf_list"),
        )
        .arg(
            Arg::with_name("vcf_list")
                .short("l")
                .long("vcf-list")
                .value_name("VCF_LIST_FILE")
                .help("File containing list of input VCF files, one per line")
                .conflicts_with("vcf"),
        )
        .arg(
            Arg::with_name("reference")
                .short("r")
                .long("reference")
                .value_name("FASTA_FILE")
                .help("Reference FASTA file")
                .required(true),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .value_name("OUTPUT_FILE")
                .help("Output file name")
                .required(true),
        )
        .arg(
            Arg::with_name("workers")
                .short("w")
                .long("workers")
                .value_name("NUM_WORKERS")
                .help("Number of worker threads")
                .required(false),
        )
        .get_matches();

    // Get the VCF files from command-line or VCF list file
    let vcf_files: Vec<PathBuf> = if let Some(vcf_list_file) = matches.value_of("vcf_list") {
        match read_vcf_list(vcf_list_file) {
            Ok(files) => files,
            Err(e) => {
                eprintln!("Error reading VCF list file {}: {}", vcf_list_file, e);
                std::process::exit(1);
            }
        }
    } else {
        matches
            .values_of("vcf")
            .unwrap_or_default()
            .map(PathBuf::from)
            .collect()
    };

    if vcf_files.is_empty() {
        eprintln!("No VCF files provided. Use -v to specify VCF files or -l to provide a VCF list file.");
        std::process::exit(1);
    }

    let reference_file = matches.value_of("reference").unwrap();
    let output_file = matches.value_of("output").unwrap();

    // Get the number of worker threads
    let num_workers = matches
        .value_of("workers")
        .and_then(|w| w.parse::<usize>().ok())
        .unwrap_or_else(num_cpus::get); // Default to the number of logical CPUs

    // Initialize the Rayon thread pool with the specified number of workers
    ThreadPoolBuilder::new()
        .num_threads(num_workers)
        .build_global()
        .unwrap_or_else(|e| {
            eprintln!("Failed to build Rayon thread pool: {}", e);
            std::process::exit(1);
        });

    // Read and filter VCF files
    let variants = read_and_concatenate_vcfs(&vcf_files);

    // Extract ALT keys (unique ALT alleles)
    let alt_keys = get_alt_keys(&variants);

    // Read the reference sequence
    let reference = match read_reference_sequence(reference_file) {
        Ok(seq) => seq,
        Err(e) => {
            eprintln!("Error reading reference file {}: {}", reference_file, e);
            std::process::exit(1);
        }
    };

    // Group variants by sample
    let sample_variants_map = group_variants_by_sample(&variants);

    // Get the list of sample names
    let cov_list = get_cov_list(&variants);

    if cov_list.is_empty() {
        eprintln!("No samples found in the provided VCF files.");
        std::process::exit(1);
    }

    // Generate sample pairs (only upper triangle)
    let sample_pairs = generate_sample_pairs(&cov_list);

    // Set up the progress bar
    let pb = ProgressBar::new(sample_pairs.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({percent}%) {msg}")
            .progress_chars("#>-"),
    );

    // Compute pairwise distances in parallel and collect the results
    let pairwise_results: Vec<((usize, usize), f64)> = sample_pairs
        .par_iter()
        .map(|&(i, j)| {
            let sample1 = &cov_list[i];
            let sample2 = &cov_list[j];

            let sample_variants1 = match sample_variants_map.get(sample1) {
                Some(v) => v,
                None => {
                    eprintln!("Sample {} not found in variants map.", sample1);
                    &Vec::new()
                }
            };
            let sample_variants2 = match sample_variants_map.get(sample2) {
                Some(v) => v,
                None => {
                    eprintln!("Sample {} not found in variants map.", sample2);
                    &Vec::new()
                }
            };

            let sum_fst = get_dif_n(sample_variants1, sample_variants2, &reference, &alt_keys);

            pb.inc(1);
            ((i, j), sum_fst)
        })
        .collect();

    pb.finish_with_message("Calculation complete");

    // Initialize the distance matrix with zeros on the diagonal
    let num_samples = cov_list.len();
    let mut distance_matrix: Vec<Vec<f64>> = vec![vec![0.0; num_samples]; num_samples];

    // Fill the upper triangle with computed Fst values
    for ((i, j), fst) in pairwise_results {
        distance_matrix[i][j] = fst;
    }

    // Fill the lower triangle by mirroring the upper triangle
    for i in 0..num_samples {
        for j in 0..i {
            distance_matrix[i][j] = distance_matrix[j][i];
        }
    }

    // Write the distance matrix to the output file
    if let Err(e) = write_distance_matrix(output_file, &distance_matrix, &cov_list) {
        eprintln!("Error writing output file {}: {}", output_file, e);
        std::process::exit(1);
    }

    println!("Distance matrix successfully written to {}", output_file);
}

/// Reads a file containing a list of VCF file paths, one per line.
/// Ignores empty lines and lines starting with '#' (comments).
fn read_vcf_list(vcf_list_file: &str) -> Result<Vec<PathBuf>, std::io::Error> {
    let file = File::open(vcf_list_file)?;
    let reader = BufReader::new(file);
    let mut vcf_files = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        vcf_files.push(PathBuf::from(trimmed));
    }

    Ok(vcf_files)
}

/// Parses a VCF file and returns a vector of Variants.
fn parse_vcf(file: &PathBuf) -> Vec<Variant> {
    let sample_id = file
        .file_name()
        .and_then(|s| s.to_str())
        .map(|s| s.split('.').next().unwrap_or(""))
        .unwrap_or("")
        .to_string();

    let f = match File::open(file) {
        Ok(file) => file,
        Err(e) => {
            eprintln!("Error opening VCF file {}: {}", file.display(), e);
            return Vec::new();
        }
    };
    let reader = BufReader::new(f);

    let mut data = Vec::new();

    for line in reader.lines() {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                eprintln!("Error reading line in {}: {}", file.display(), e);
                continue;
            }
        };

        if line.starts_with("##") {
            continue;
        }

        if line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            continue; // Skip incomplete lines
        }

        let pos = match fields[1].parse::<usize>() {
            Ok(p) => p,
            Err(_) => continue,
        };
        let ref_allele = fields[3];
        let alt_allele = fields[4];
        let format_fields = fields[8];
        let sample_data = fields[9];

        // Parse FORMAT and sample data
        let format_keys: Vec<&str> = format_fields.split(':').collect();
        let format_values: Vec<&str> = sample_data.split(':').collect();

        let format_dict: HashMap<&str, &str> = format_keys
            .iter()
            .zip(format_values.iter())
            .map(|(k, v)| (*k, *v))
            .collect();

        // Extract fields we need
        let total_dp = format_dict
            .get("DP")
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(0);
        let ref_dp = format_dict
            .get("RD")
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(0);
        let alt_dp = format_dict
            .get("AD")
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(0);
        let alt_freq_str = format_dict.get("FREQ").unwrap_or(&"0%");
        let alt_freq = alt_freq_str.trim_end_matches('%').parse::<f64>().unwrap_or(0.0) / 100.0;
        let alt_rv = format_dict
            .get("ADR")
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(0);

        let variant = Variant {
            region: sample_id.clone(),
            pos,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            total_dp,
            ref_dp,
            alt_dp,
            alt_freq,
            alt_rv,
        };

        data.push(variant);
    }

    data
}

/// Reads and filters VCF files, returning a concatenated list of Variants.
fn read_and_concatenate_vcfs(file_list: &[PathBuf]) -> Vec<Variant> {
    file_list
        .par_iter()
        .flat_map(|file| read_and_filter_vcf(file))
        .collect()
}

/// Reads and filters a single VCF file, returning a vector of Variants.
fn read_and_filter_vcf(file: &PathBuf) -> Vec<Variant> {
    let variants = parse_vcf(file);

    let filtered_variants: Vec<Variant> = variants
        .into_iter()
        .filter(|v| {
            let depth_filter = v.total_dp >= 30;
            let freq_filter = v.alt_freq >= 0.05;
            let num_reads = v.alt_dp >= 2;
            let num_reads_rev = v.alt_rv >= 2;
            depth_filter && freq_filter && num_reads && num_reads_rev
        })
        .collect();

    filtered_variants
}

/// Extracts unique ALT alleles from the list of Variants.
fn get_alt_keys(variants: &[Variant]) -> Vec<String> {
    let mut alt_set: HashSet<String> = HashSet::new();
    for variant in variants {
        alt_set.insert(variant.alt_allele.clone());
    }
    let mut alt_keys: Vec<String> = alt_set.into_iter().collect();
    alt_keys.sort();
    alt_keys
}

/// Reads the reference sequence from a FASTA file.
fn read_reference_sequence(reference_file: &str) -> Result<Vec<char>, std::io::Error> {
    let f = File::open(reference_file)?;
    let reader = BufReader::new(f);
    let mut sequence = Vec::new();
    let mut seq_started = false;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            seq_started = true;
            continue;
        }
        if seq_started {
            sequence.extend(line.trim().chars());
        }
    }

    Ok(sequence)
}

/// Groups variants by sample.
fn group_variants_by_sample(variants: &[Variant]) -> HashMap<String, Vec<Variant>> {
    let mut sample_variants_map: HashMap<String, Vec<Variant>> = HashMap::new();
    for variant in variants {
        sample_variants_map
            .entry(variant.region.clone())
            .or_insert_with(Vec::new)
            .push(variant.clone());
    }
    sample_variants_map
}

/// Retrieves a sorted list of unique sample names from the variants.
fn get_cov_list(variants: &[Variant]) -> Vec<String> {
    let mut sample_set = HashSet::new();
    for variant in variants {
        sample_set.insert(variant.region.clone());
    }
    let mut cov_list: Vec<String> = sample_set.into_iter().collect();
    cov_list.sort();
    cov_list
}

/// Generates all unique sample pairs for upper triangle calculation.
fn generate_sample_pairs(samples: &[String]) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for i in 0..samples.len() {
        for j in (i + 1)..samples.len() {
            pairs.push((i, j));
        }
    }
    pairs
}

/// Retrieves allele frequencies for a given position.
fn get_pos_tup(
    sample_variants: &[Variant],
    pos: usize,
    reference: &[char],
    alt_keys: &[String],
) -> Vec<f64> {
    let mut freq: HashMap<String, f64> = alt_keys.iter().map(|key| (key.clone(), 0.0)).collect();

    if let Some(variant) = sample_variants.iter().find(|v| v.pos == pos + 1) {
        let alt_allele = &variant.alt_allele;
        if freq.contains_key(alt_allele) {
            freq.insert(alt_allele.clone(), variant.alt_freq);
        }
    }

    let sum_freq: f64 = freq.values().sum();
    let ref_base = reference[pos].to_string();
    let ref_freq = 1.0 - sum_freq;
    *freq.entry(ref_base).or_insert(0.0) += ref_freq;

    alt_keys
        .iter()
        .map(|key| *freq.get(key).unwrap_or(&0.0))
        .collect()
}

/// Calculates heterozygosity based on allele frequencies.
fn heterozygosity(freqs: &[f64]) -> f64 {
    1.0 - freqs.iter().map(|f| f * f).sum::<f64>()
}

/// Calculates hs and ht for two samples at a given position.
fn calc_heterozygosities(
    sample_variants1: &[Variant],
    sample_variants2: &[Variant],
    pos: usize,
    reference: &[char],
    alt_keys: &[String],
) -> (f64, f64) {
    let freqs1 = get_pos_tup(sample_variants1, pos, reference, alt_keys);
    let freqs2 = get_pos_tup(sample_variants2, pos, reference, alt_keys);

    let hs1 = heterozygosity(&freqs1);
    let hs2 = heterozygosity(&freqs2);
    let hs = (hs1 + hs2) / 2.0;

    let total_freqs: Vec<f64> = freqs1
        .iter()
        .zip(freqs2.iter())
        .map(|(f1, f2)| (f1 + f2) / 2.0)
        .collect();

    let ht = heterozygosity(&total_freqs);

    (hs, ht)
}

/// Calculates Fst using Weir and Cockerham's method.
fn calc_fst_weir_cockerham(hs: f64, ht: f64) -> f64 {
    if ht != 0.0 {
        (ht - hs) / ht
    } else {
        0.0
    }
}

/// Calculates the sum of Fst values for a sample pair across all positions.
fn get_dif_n(
    sample_variants1: &[Variant],
    sample_variants2: &[Variant],
    reference: &[char],
    alt_keys: &[String],
) -> f64 {
    let length = reference.len();

    // Parallelize over positions
    (0..length)
        .into_par_iter()
        .map(|pos| {
            let (hs, ht) =
                calc_heterozygosities(sample_variants1, sample_variants2, pos, reference, alt_keys);
            calc_fst_weir_cockerham(hs, ht)
        })
        .sum()
}

/// Writes the distance matrix to the specified output file.
fn write_distance_matrix(
    output_file: &str,
    distance_matrix: &[Vec<f64>],
    cov_list: &[String],
) -> std::io::Result<()> {
    let mut file = File::create(output_file)?;

    // Write header
    let header = format!("\t{}", cov_list.join("\t"));
    writeln!(file, "{}", header)?;

    for (i, row) in distance_matrix.iter().enumerate() {
        let mut line = vec![cov_list[i].clone()];
        for &value in row {
            let s = format!("{:.6}", value);
            line.push(s);
        }
        writeln!(file, "{}", line.join("\t"))?;
    }

    Ok(())
}
