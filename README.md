<p align="center">
  <a href="https://github.com/PathoGenOmics-Lab/fstic">
    <img src="https://github.com/PathoGenOmics-Lab/fstic/blob/main/.github/logos/fstic.png" height="300" alt="fstic">
  </a>
</p>

# Fstic: Allele‑Frequency‑based Genetic distance calculator
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/PathoGenOmics-Lab/fstic/blob/main/LICENSE)
[![distree](https://img.shields.io/badge/fstic-rust-%23ff8000)](https://github.com/PathoGenOmics-Lab/fstic)
[![Anaconda-Server Badge](https://img.shields.io/conda/dn/bioconda/fstic.svg)](https://anaconda.org/bioconda/fstic)
[![Anaconda-Version Badge](https://anaconda.org/bioconda/fstic/badges/version.svg)](https://anaconda.org/bioconda/fstic)
[![PGO](https://img.shields.io/badge/PathoGenOmics-lab-red?)](https://github.com/PathoGenOmics-Lab)
[![DOI](https://img.shields.io/badge/doi-10.5281%2Fzenodo.16813662-%23ff0077)](https://doi.org/10.5281/zenodo.16813662)


__Paula Ruiz-Rodriguez<sup>1</sup>__ 
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. I<sup>2</sup>SysBio, University of Valencia-CSIC, FISABIO Joint Research Unit Infection and Public Health, Valencia, Spain </sub>  


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
# 1. Install your program:

# Using conda
conda install -c bioconda fstic
or
# Using mamba
mamba install -c bioconda fstic
or
# Generate your binary (requires Rust toolchain)
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
| `--min-alt-rev-reads <INT>` | 2                        | Filter out variants with fewer supporting alt reads on the reverse strand.                           |
| `--pass-only`            | *off*                       | Keep only variants whose VCF `FILTER` column is `PASS` or `.`.                                       |
| `--workers <INT>`        | *all logical cores*         | Number of threads to spawn.                                                                          |
| `--help`                 | –                           | Print the full help message.                                                                         |

---

## Guide to Distance Formulas

### Metrics for Population Differentiation

| Name (`--formula`) | Global Formula                                                                                                 | Notes & Recommended Use                                                                                                          |
| ------------------ | -------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| **GST**            | \$G\_{ST} = \dfrac{\sum\_l (H\_{T,l}-H\_{S,l})}{\sum\_l H\_{T,l}}\$ <br>*(ratio‑of‑sums)* | Classic overall differentiation (Nei 1973). Bounded in [0, 1]; unaffected by `--normalize`. |
| **FST**            | \$F\_{ST} = \sum\_l \dfrac{H\_{T,l}-H\_{S,l}}{H\_{T,l}}\$ <br>*(sum of per‑locus \$G\_{ST}\$)* | Per‑locus Nei \$G\_{ST}\$ summed over sites, so it grows with the number of loci; pass `--normalize` for the mean. Not the Weir & Cockerham \$\theta\$ estimator. Default. |
| **Jost’s D**       | \$D = \dfrac{n}{n-1} \cdot \dfrac{H\_T - H\_S}{1 - H\_S}\$                                                                          | Measures the fraction of allelic diversity that is partitioned among populations; less sensitive to within‑population variation. |

### Metrics for Phylogenetic / Divergence Analysis

| Name (`--formula`)         | Formula                                                       | Notes                                                                         |
| -------------------------- | ------------------------------------------------------------- | ----------------------------------------------------------------------------- |
| **Reynolds**               | \$D\_R = -\ln(1 - \theta)\$, with <br>\$\theta = \dfrac{\sum\_l \sum\_i (p\_i-q\_i)^2 / 2}{\sum\_l \left(1 - \sum\_i p\_i q\_i\right)}\$ | Reynolds, Weir & Cockerham (1983) coancestry coefficient. Linear with drift time for recently diverged populations. |
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
* SNPs only. Indels, multi‑allelic records (`ALT` with a comma), `ALT=.` and `ALT=*` are skipped, and the counts are reported at the end of the run.
* `--reference` **not** required in VCF mode: REF allele is part of the file.

**Example 1 — VCF with all filtering fields**

```vcf
##fileformat=VCFv4.2
#CHROM POS  ID REF ALT QUAL FILTER INFO FORMAT        sample_A
chr1    100 .  A   T   .    .      .    GT:DP:AD:FREQ 0/1:50:24:48.0%
chr1    250 .  C   G   .    .      .    GT:DP:AD:FREQ 0/1:45:10:22.2%
```

`AD` is read as the standard `ref,alt` pair when it carries a comma, and as a bare
alt count otherwise, which is what VarScan writes. When `FREQ` is absent it is
computed from `AD/DP`.

**Example 2 — VCF missing DP/AD but with FREQ only**

```vcf
##fileformat=VCFv4.2
#CHROM POS  ID REF ALT QUAL FILTER INFO FORMAT sample_B
chr1    100 .  A   T   .    .      .    GT:FREQ 0/1:85.0%
```

---

### Table (`.csv`, `.tsv`, `.tab`)

* **Required columns:** `sample`, `position`, `sequence` (alt allele), `frequency`.
* **Recommended:** `ref_allele`. Without it the reference base is taken from `--reference`.
* **Optional filtering columns:** `total_dp`, `alt_dp`, `alt_rv`.
* Optional `chrom` column. Either every input table has it or none does; mixing the two is rejected, since the same site would otherwise be counted twice.
* Column names are case‑insensitive; delimiter auto‑detected from extension.
* `frequency` accepts a proportion (`0.5`) or a percentage (`50%`), as in VCF `FREQ`. A bare number is always a proportion.

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
<h2 id="contributors" align="center">

✨ [Contributors]((https://github.com/PathoGenOmics-Lab/AMAP/graphs/contributors))
</h2>

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<div align="center">
fstic is developed with ❤️ by:
<table>
  <tr>
    <td align="center">
      <a href="https://github.com/paururo">
        <img src="https://avatars.githubusercontent.com/u/50167687?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Paula Ruiz-Rodriguez</b></sub>
      </a>
      <br />
      <a href="" title="Code">💻</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Data">🔣</a>
      <a href="" title="Desing">🎨</a>
      <a href="" title="Tool">🔧</a>
    </td> 
    <td align="center">
      <a href="https://github.com/mireiacoscolla">
        <img src="https://avatars.githubusercontent.com/u/29301737?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Mireia Coscolla</b></sub>
      </a>
      <br />
      <a href="https://www.uv.es/instituto-biologia-integrativa-sistemas-i2sysbio/es/investigacion/proyectos/proyectos-actuales/mol-tb-host-1286169137294/ProjecteInves.html?id=1286289780236" title="Funding/Grant Finders">🔍</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Mentoring">🧑‍🏫</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="User Testing">📓</a>
    </td> 
  </tr>
</table>

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification ([emoji key](https://allcontributors.org/docs/en/emoji-key)).

<!-- markdownlint-restore -->

---



