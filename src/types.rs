use std::collections::HashMap;

/// Genomic position: chromosome + 1-based coordinate.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct GenomicPos {
    pub chrom: String,
    pub pos: usize,
}

#[derive(Clone, Debug, Default)]
pub struct SiteData {
    pub reference_allele: String,
    pub freqs: AlleleFrequencies,
}

pub type AlleleFrequencies = HashMap<String, f64>;
pub type PositionalData = HashMap<GenomicPos, SiteData>;
pub type SampleVariants = HashMap<String, PositionalData>;

pub struct VcfVariant {
    pub sample: String,
    pub chrom: String,
    pub pos: usize,
    pub ref_allele: String,
    pub alt_allele: String,
    pub total_dp: Option<u32>,
    pub alt_dp: Option<u32>,
    pub alt_freq: f64,
    pub alt_rv: Option<u32>,
}

#[derive(Debug, Clone, Copy)]
pub struct FilterCriteria {
    pub min_depth: u32,
    pub min_freq: f64,
    pub min_alt_reads: u32,
    pub min_alt_rev_reads: u32,
    pub pass_only: bool,
}

#[derive(Debug, serde::Deserialize)]
pub struct TableInputRow {
    pub sample: String,
    #[serde(default = "default_chrom")]
    pub chrom: String,
    pub position: usize,
    pub sequence: String,
    pub frequency: f64,
    pub ref_allele: Option<String>,
    pub total_dp: Option<u32>,
    pub alt_dp: Option<u32>,
    pub alt_rv: Option<u32>,
}

fn default_chrom() -> String {
    ".".to_string()
}
