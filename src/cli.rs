use anyhow::{bail, Context, Result};
use clap::Parser;
use rayon::ThreadPoolBuilder;
use std::path::PathBuf;

#[derive(Debug, Clone, Copy)]
pub enum InputMode {
    Vcf,
    Table,
}

#[derive(Parser)]
#[command(
    name = "fstic",
    version,
    about = "Calculates pairwise genetic distances between samples using standard estimators.",
    author = "Paula Ruiz-Rodriguez <paula.ruiz.rodriguez@csic.es>"
)]
pub struct Cli {
    /// One or more input VCF files
    #[arg(long = "vcf", value_name = "VCF_FILES", num_args = 1.., conflicts_with_all = ["vcf_list", "table", "table_list"])]
    pub vcf: Option<Vec<PathBuf>>,

    /// A file containing a list of input VCF files
    #[arg(long = "vcf-list", value_name = "VCF_LIST", conflicts_with_all = ["vcf", "table", "table_list"])]
    pub vcf_list: Option<PathBuf>,

    /// One or more input table files (format detected by extension: .csv, .tsv, .tab)
    #[arg(long = "table", value_name = "TABLE_FILES", num_args = 1.., conflicts_with_all = ["vcf", "vcf_list", "table_list"])]
    pub table: Option<Vec<PathBuf>>,

    /// A file containing a list of input table files
    #[arg(long = "table-list", value_name = "TABLE_LIST", conflicts_with_all = ["vcf", "vcf_list", "table"])]
    pub table_list: Option<PathBuf>,

    /// Reference FASTA file (required for table mode, optional for VCF)
    #[arg(short = 'r', long = "reference", value_name = "FASTA_FILE")]
    pub reference: Option<PathBuf>,

    /// Output file name for the distance matrix
    #[arg(short = 'o', long = "output", value_name = "OUTPUT_FILE")]
    pub output: PathBuf,

    /// The distance formula to use
    #[arg(short = 'f', long = "formula", value_name = "FORMULA", default_value = "fst",
          value_parser = ["fst", "gst", "nei", "chord", "bray-curtis", "jost_d", "reynolds", "rogers"])]
    pub formula: String,

    /// Normalize by the number of loci (affects: fst, chord, bray-curtis, jost_d)
    #[arg(long = "normalize")]
    pub normalize: bool,

    /// Minimum total read depth to keep a variant
    #[arg(long = "min-depth", value_name = "DP", default_value = "30")]
    pub min_depth: u32,

    /// Minimum alternate allele frequency to keep a variant
    #[arg(long = "min-af", value_name = "FREQ", default_value = "0.05")]
    pub min_af: f64,

    /// Minimum number of alternate allele reads to keep a variant
    #[arg(long = "min-alt-reads", value_name = "AD", default_value = "2")]
    pub min_alt_reads: u32,

    /// Minimum number of alternate allele reverse reads to keep a variant
    #[arg(long = "min-alt-rev-reads", value_name = "ADR", default_value = "2")]
    pub min_alt_rev_reads: u32,

    /// Number of worker threads (default: all available cores)
    #[arg(short = 'w', long = "workers", value_name = "NUM_WORKERS")]
    pub workers: Option<usize>,
}

impl Cli {
    pub fn get_input_files(&self) -> Result<(InputMode, Vec<PathBuf>)> {
        if let Some(ref files) = self.vcf {
            return Ok((InputMode::Vcf, files.clone()));
        }
        if let Some(ref list_path) = self.vcf_list {
            let path_str = list_path.to_str().unwrap_or("");
            let files = super::io::read_file_list(path_str)
                .with_context(|| format!("Failed to read VCF list file: {}", path_str))?;
            return Ok((InputMode::Vcf, files));
        }
        if let Some(ref files) = self.table {
            return Ok((InputMode::Table, files.clone()));
        }
        if let Some(ref list_path) = self.table_list {
            let path_str = list_path.to_str().unwrap_or("");
            let files = super::io::read_file_list(path_str)
                .with_context(|| format!("Failed to read table list file: {}", path_str))?;
            return Ok((InputMode::Table, files));
        }
        bail!("No input files provided. Use --vcf, --vcf-list, --table, or --table-list.");
    }

    pub fn configure_thread_pool(&self) -> Result<()> {
        let num_workers = self.workers.unwrap_or_else(num_cpus::get);
        ThreadPoolBuilder::new()
            .num_threads(num_workers)
            .build_global()
            .context("Failed to initialize thread pool")?;
        Ok(())
    }
}
