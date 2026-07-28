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
    /// `None` when the table carries no `chrom` column. Kept distinct from a real
    /// contig name so that mixing tables with and without one can be detected.
    #[serde(default)]
    pub chrom: Option<String>,
    pub position: usize,
    pub sequence: String,
    #[serde(deserialize_with = "de_frequency")]
    pub frequency: f64,
    pub ref_allele: Option<String>,
    pub total_dp: Option<u32>,
    pub alt_dp: Option<u32>,
    pub alt_rv: Option<u32>,
}

/// Accepts a proportion (`0.125`) or a percentage (`12.5%`), matching the VCF FREQ
/// field. A bare number is always a proportion: guessing from magnitude would turn a
/// column mistakenly holding read counts into plausible-looking frequencies.
fn de_frequency<'de, D>(deserializer: D) -> Result<f64, D::Error>
where
    D: serde::Deserializer<'de>,
{
    use serde::Deserialize;
    let raw = String::deserialize(deserializer)?;
    let s = raw.trim();
    let parsed = match s.strip_suffix('%') {
        Some(stripped) => stripped.trim().parse::<f64>().map(|f| f / 100.0),
        None => s.parse::<f64>(),
    };
    parsed.map_err(|_| serde::de::Error::custom(format!("invalid frequency '{}'", raw)))
}
