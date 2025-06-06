use std::collections::HashMap;

#[derive(Clone, Debug, Default)]
pub struct SiteData {
    pub reference_allele: String,
    pub freqs: AlleleFrequencies,
}

pub type Position = usize;
pub type Allele = String;
pub type Frequency = f64;
pub type AlleleFrequencies = HashMap<Allele, Frequency>;
pub type PositionalData = HashMap<Position, SiteData>;
pub type SampleVariants = HashMap<String, PositionalData>;

pub enum Calculation {
    Fst,
    Gst,
    Nei,
    Chord,
    BrayCurtis,
    JostD,
    Reynolds,
    Rogers,
}

pub struct VcfVariant {
    pub sample: String,
    pub pos: usize,
    pub ref_allele: String,
    pub alt_allele: String,
    pub total_dp: Option<u32>,
    pub alt_dp: Option<u32>,
    pub alt_freq: Option<f64>,
    pub alt_rv: Option<u32>,
}

#[derive(Debug, Clone, Copy)]
pub struct FilterCriteria {
    pub min_depth: u32,
    pub min_freq: f64,
    pub min_alt_reads: u32,
    pub min_alt_rev_reads: u32,
}

#[derive(Debug, serde::Deserialize)]
pub struct TableInputRow {
    pub sample: String,
    pub position: usize,
    pub sequence: String,
    pub frequency: f64,
    pub ref_allele: Option<String>,
    pub total_dp: Option<u32>,
    pub alt_dp: Option<u32>,
    pub alt_rv: Option<u32>,
}
