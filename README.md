<p align="center">
  <a href="https://github.com/PathoGenOmics-Lab/fstigo">
    <img src="https://github.com/PathoGenOmics-Lab/fstigo/blob/main/.github/logs/fstigo.png" height="250" alt="fstic">
  </a>
</p>

# Fstic: Allele‑Frequency‑based Genetic distance calculator

---

## Overview

Fstic is a high‑performance command‑line tool written in Rust that calculates pairwise genetic distances from variant data. It ingests single‑sample VCF files **or** allele‑frequency tables and outputs an *N × N* distance matrix that summarises genetic differentiation among samples. The code is fully parallelised and scales to whole‑genome datasets on modern multi‑core CPUs.

---

## Key Features

* **Multiple Distance Metrics** · Eight standard estimators (FST, GST, Jost’s *D*, Reynolds, Nei, Cavalli‑Sforza chord, Rogers, Bray‑Curtis) let you view your data from complementary theoretical angles.
* **Flexible Input** · Accepts raw VCFs or pre‑processed tables (`.csv`, `.tsv`, `.tab`). You can pass file *lists* for convenience.
* **Smart Filtering** · Configurable filters on depth, allele frequency and allele counts for both VCF and table inputs ensure data quality.
* **Configurable Output** · Add `--normalize` to divide cumulative distances by the number of loci (where applicable).
* **Optimised for Speed** · Work is automatically distributed across all logical CPU cores (override with `--workers`).
* **Transparent Reporting** · Live progress bar, ETA and a complete log of the filters applied.

---

## Quick Start

```bash
# 1. Build (requires Rust toolchain)
git clone https://github.com/<your-org>/fstic.git
cd fstic
cargo build --release  # binary at ./target/release/fstic

# 2. Run examples

# Example A: basic FST from three VCFs
./fstic \
  --vcf sample1.vcf sample2.vcf sample3.vcf \
  --output distances.csv

# Example B: Bray–Curtis from a list of tables, eight threads
./fstic \
  --table-list path_to_tables.txt \
  --output bray_curtis.csv \
  --formula bray-curtis \
  --workers 8

# Example C: stringent filters + normalised FST
./fstic \
  --table *.tsv \
  --reference ref.fa \
  --output fst_norm.tsv \
  --formula fst \
  --normalize \
  --min-depth 50 \
  --min-af 0.01 \
  --min-alt-reads 5
```

---

## Command‑Line Reference

| Flag                     | Default                     | Description                                                                                          |
| ------------------------ | --------------------------- | ---------------------------------------------------------------------------------------------------- |
| `--vcf / --vcf-list`     | –                           | One or more VCF files or a file containing paths to them.                                            |
| `--table / --table-list` | –                           | One or more table files or a file with paths.                                                        |
| `--reference <FASTA>`    | *required for table inputs* | Reference genome used to infer missing reference alleles. Not needed for VCF mode.                   |
| `--output <FILE>`        | *required*                  | Destination file for the distance matrix.                                                            |
| `--formula <METRIC>`     | `fst`                       | Distance metric to compute (see below).                                                              |
| `--normalize`            | *off*                       | Divide cumulative distances by the number of loci (affects `fst`, `bray-curtis`, `chord`, `jost_d`). |
| `--min-depth <INT>`      | 30                          | Filter out variants whose total depth < this value.                                                  |
| `--min-af <FLOAT>`       | 0.05                        | Filter out variants whose alt‑allele frequency < this value.                                         |
| `--min-alt-reads <INT>`  | 2                           | Filter out variants with fewer supporting alt reads.                                                 |
| `--workers <INT>`        | *all logical cores*         | Number of threads to spawn.                                                                          |
| `--help`                 | –                           | Print the full help message.                                                                         |

---

## Guide to Distance Formulas

### Metrics for Population Differentiation

| Name (`--formula`) | Global Formula                                                                                                 | Notes & Recommended Use                                                                                                          |
| ------------------ | -------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| **GST**            | \$G\_{ST} = \dfrac{H\_T - H\_S}{H\_T}\$                                                                        | Classic overall differentiation (Nei 1973).                                                                                      |
| **FST**            | \$\theta = \dfrac{\sum\_l (H\_{T,l}-H\_{S,l})}{\sum\_l H\_{T,l}}\$ <br>*(Weir & Cockerham 1984 ratio‑of‑sums)* | Default for relative differentiation; per‑locus estimates are also available.                                                    |
| **Jost’s D**       | \$D = \dfrac{H\_T - H\_S}{1 - H\_S}\$                                                                          | Measures the fraction of allelic diversity that is partitioned among populations; less sensitive to within‑population variation. |

### Metrics for Phylogenetic / Divergence Analysis

| Name (`--formula`)         | Formula                                                       | Notes                                                                         |
| -------------------------- | ------------------------------------------------------------- | ----------------------------------------------------------------------------- |
| **Reynolds**               | \$D\_R = -\ln(1 - \theta)\$                                   | Linear with drift time for recently diverged populations.                     |
| **Nei’s D**                | \$D = -\ln \left( \dfrac{J\_{xy}}{\sqrt{J\_x J\_y}} \right)\$ | Effective over long time‑scales; \$J\$ = probability of allele identity.      |
| **Cavalli‑Sforza “Chord”** | \$D\_{CH} = \sqrt{,2\bigl(1-\sum\_i \sqrt{p\_i q\_i}\bigr)}\$ | Geometric distance satisfying triangle inequality; useful for tree inference. |
| **Rogers**                 | \$D\_{R} = \sqrt{\dfrac{\sum\_i (p\_i-q\_i)^2}{2L}}\$         | Euclidean‑based distance bounded between 0 and 1.                             |

### Metric for Allele‑Frequency Profile Comparison

| Name (`--formula`) | Formula                                                 | Notes                                                         |
| ------------------ | ------------------------------------------------------- | ------------------------------------------------------------- |
| **Bray‑Curtis**    | \$BC = \tfrac{1}{2} \sum\_i \lvert p\_i - q\_i \rvert\$ | **Dissimilarity** (does not satisfy the triangle inequality). |

*Variables:* \$p\_i\$, \$q\_i\$ = allele frequencies in populations *x*, *y*; \$H\_S\$ and \$H\_T\$ = within‑population and total heterozygosity; \$L\$ = number of loci.

---

## Input Formats

### VCF

* **One file per sample** (single‑sample VCF). Sample name is taken from the filename.
* Requires `FORMAT/FREQ` field (supports decimal **or** percentage). Optionally uses `DP`, `AD`, `ADR` for filtering.
* Indels fully supported.
* `--reference` **not** required in VCF mode: REF allele is part of the file.

**Example 1 — VCF with all filtering fields**

```vcf
##fileformat=VCFv4.2
#CHROM POS  ID REF ALT QUAL FILTER INFO FORMAT        sample_A
chr1    100 .  A   T   .    .      .    GT:DP:AD:FREQ 0/1:50:24:48.0%
chr1    250 .  C   CAA .    .      .    GT:DP:AD:FREQ 0/1:45:10:22.2%
```

**Example 2 — VCF missing DP/AD but with FREQ only**

```vcf
##fileformat=VCFv4.2
#CHROM POS  ID REF ALT QUAL FILTER INFO FORMAT sample_B
chr1    100 .  A   T   .    .      .    GT:FREQ 0/1:85.0%
```

---

### Table (`.csv`, `.tsv`, `.tab`)

* **Required columns:** `sample`, `position`, `sequence` (alt allele), `frequency`.
* **Recommended:** `ref_allele` (crucial for indels).
* **Optional filtering columns:** `total_dp`, `alt_dp`, `alt_rv`.
* Column names are case‑insensitive; delimiter auto‑detected from extension.

**Example 1 — CSV with all columns**

```csv
sample,position,ref_allele,sequence,frequency,total_dp,alt_dp
sample_A,100,A,T,0.5,50,24
sample_B,100,A,T,0.8,60,48
```

**Example 2 — TSV with only required columns**

```tsv
sample	position	sequence	frequency
sample_C	550	T	0.12
sample_D	550	A	0.95
```

---

## Minor Clarifications & Defaults

* **Threads:** By default Fstic launches one worker per *logical* CPU core. Override with `--workers N`.
* **Frequencies:** The `frequency` field can be a proportion (0.125) **or** a percentage (12.5 %). Both are auto‑detected.
* **Normalisation:** Adding `--normalize` turns cumulative distances into per‑locus means for metrics where that is meaningful.

---



