use clap::Parser;
use rayon::ThreadPoolBuilder;
use std::path::PathBuf;

#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum Formula {
    Fst,
    Gst,
    Nei,
    Chord,
    #[value(name = "bray-curtis")]
    BrayCurtis,
    #[value(name = "jost_d")]
    JostD,
    Reynolds,
    Rogers,
}

pub enum InputMode {
    Vcf,
    Table,
}

/// High-performance pairwise genetic distance calculator from allele-frequency data.
#[derive(Parser, Debug)]
#[command(name = "fstic", version, about, author = "Paula Ruiz-Rodriguez <paula.ruiz.rodriguez@csic.es>")]
pub struct Args {
    /// One or more input VCF files.
    #[arg(long = "vcf", num_args = 1.., value_name = "VCF_FILES", group = "input")]
    pub vcf: Option<Vec<PathBuf>>,

    /// A file containing a list of input VCF paths.
    #[arg(long = "vcf-list", value_name = "VCF_LIST", group = "input")]
    pub vcf_list: Option<PathBuf>,

    /// One or more input table files (.csv, .tsv, .tab).
    #[arg(long = "table", num_args = 1.., value_name = "TABLE_FILES", group = "input")]
    pub table: Option<Vec<PathBuf>>,

    /// A file containing a list of input table paths.
    #[arg(long = "table-list", value_name = "TABLE_LIST", group = "input")]
    pub table_list: Option<PathBuf>,

    /// Reference FASTA file (required for table inputs; optional for VCF).
    #[arg(short = 'r', long = "reference", value_name = "FASTA_FILE")]
    pub reference: Option<PathBuf>,

    /// Output file name for the distance matrix.
    #[arg(short = 'o', long = "output", value_name = "OUTPUT_FILE")]
    pub output: PathBuf,

    /// Distance formula to use.
    #[arg(short = 'f', long = "formula", default_value = "fst")]
    pub formula: Formula,

    /// Normalize by the number of loci (affects: fst, chord, bray-curtis, jost_d).
    #[arg(long = "normalize")]
    pub normalize: bool,

    /// Minimum total read depth to keep a variant.
    #[arg(long = "min-depth", default_value_t = 30)]
    pub min_depth: u32,

    /// Minimum alternate allele frequency to keep a variant.
    #[arg(long = "min-af", default_value_t = 0.05)]
    pub min_af: f64,

    /// Minimum number of alternate allele reads to keep a variant.
    #[arg(long = "min-alt-reads", default_value_t = 2)]
    pub min_alt_reads: u32,

    /// Minimum number of alternate allele reverse reads to keep a variant.
    #[arg(long = "min-alt-rev-reads", default_value_t = 2)]
    pub min_alt_rev_reads: u32,

    /// Only keep variants that PASS all filters (VCF FILTER column).
    #[arg(long = "pass-only")]
    pub pass_only: bool,

    /// Number of worker threads (default: all available cores).
    #[arg(short = 'w', long = "workers")]
    pub workers: Option<usize>,
}

impl Args {
    /// Resolve input mode and file list.
    pub fn get_input_files(&self) -> Result<(InputMode, Vec<PathBuf>), String> {
        if self.vcf.is_some() || self.vcf_list.is_some() {
            let files = self.resolve_files(&self.vcf, &self.vcf_list)?;
            Ok((InputMode::Vcf, files))
        } else if self.table.is_some() || self.table_list.is_some() {
            if self.reference.is_none() {
                return Err(
                    "Error: --reference is required when using table input (needed to infer reference alleles).".to_string(),
                );
            }
            let files = self.resolve_files(&self.table, &self.table_list)?;
            Ok((InputMode::Table, files))
        } else {
            Err("No input files provided. Use --vcf, --vcf-list, --table, or --table-list.".to_string())
        }
    }

    fn resolve_files(
        &self,
        direct: &Option<Vec<PathBuf>>,
        list: &Option<PathBuf>,
    ) -> Result<Vec<PathBuf>, String> {
        if let Some(list_path) = list {
            crate::io::read_file_list(
                list_path
                    .to_str()
                    .ok_or_else(|| "Invalid list file path".to_string())?,
            )
            .map_err(|e| format!("Error reading list file: {}", e))
        } else if let Some(files) = direct {
            Ok(files.clone())
        } else {
            Err("No input files provided.".to_string())
        }
    }

    /// Configure the global rayon thread pool.
    pub fn configure_thread_pool(&self) {
        let num_workers = self.workers.unwrap_or_else(num_cpus::get);
        // Silently ignore if already initialized (e.g., in tests)
        let _ = ThreadPoolBuilder::new()
            .num_threads(num_workers)
            .build_global();
    }
}
