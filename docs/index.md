---
hide:
  - navigation
---

<div class="hero" markdown>
![Fstic](assets/fstic_wordmark.svg){ .hero-logo }

# Fstic { .visually-hidden }

Pairwise genetic distance matrices from allele-frequency data. Point it at a
directory of single-sample VCFs, pick an estimator, get an *N* x *N* matrix.
</div>

## What it does

Fstic reads variant calls with allele frequencies, one file per sample, and
computes the pairwise distance between every pair of samples under one of eight
standard estimators. It is built for within-host and population data where a
sample is not a single genotype but a frequency spectrum: mixed infections,
deep-sequenced isolates, pooled populations.

```bash
fstic --vcf sample*.vcf --output distances.csv --formula fst
```

<div class="grid cards" markdown>

-   :material-function-variant: __Eight estimators__

    ---

    FST, GST, Jost's *D*, Reynolds, Nei's *D*, Cavalli-Sforza chord, Rogers and
    Bray-Curtis, each with its formula written out.

    [:octicons-arrow-right-24: Choosing one](guide/metrics.md)

-   :material-file-document-multiple: __Two input formats__

    ---

    Single-sample VCFs, or frequency tables in CSV/TSV. Pass files directly or
    a file of paths.

    [:octicons-arrow-right-24: Input formats](getting-started/inputs.md)

-   :material-filter-variant: __Filtering that reports itself__

    ---

    Depth, frequency and strand-support thresholds, with a count of everything
    dropped and why.

    [:octicons-arrow-right-24: Filtering](guide/filtering.md)

-   :material-speedometer: __Parallel and reproducible__

    ---

    Work spread across all cores, with output that does not change when the
    core count does.

    [:octicons-arrow-right-24: Reliability](how-it-works/reliability.md)

</div>

## Install

=== "conda"

    ```bash
    conda install -c bioconda fstic
    ```

=== "mamba"

    ```bash
    mamba install -c bioconda fstic
    ```

=== "from source"

    ```bash
    git clone https://github.com/PathoGenOmics-Lab/fstic.git
    cd fstic
    cargo build --release
    ```

[Getting started :octicons-arrow-right-24:](getting-started/installation.md){ .md-button .md-button--primary }
[Quick start :octicons-arrow-right-24:](getting-started/quickstart.md){ .md-button }

## What a distance means here

A sample is a vector of allele frequencies per site, not a genotype call. Where
a sample has no record at a site, it is taken to be homozygous reference there,
which is the usual convention for variant-only VCFs. The set of loci is the
union of variant sites across all inputs.

That last point has a consequence worth stating: with the cumulative estimators
a locus at which neither member of a pair carries a variant contributes nothing,
so adding an unrelated sample to a run does not move the distances between the
others. Nei's *D* and Rogers are defined over the whole locus set and do move,
by design.

[:octicons-arrow-right-24: The model in full](how-it-works/model.md)

## Citation

Fstic is developed at [I2SysBio](https://www.uv.es/instituto-biologia-integrativa-sistemas-i2sysbio/)
(University of Valencia-CSIC) by Paula Ruiz-Rodriguez and Mireia Coscolla.

[:octicons-arrow-right-24: How to cite](about/citation.md)
