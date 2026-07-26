<p align="center">
  <a href="https://pathogenomics-lab.github.io/fstic/">
    <img src="https://github.com/PathoGenOmics-Lab/fstic/blob/main/.github/logos/fstic.png" height="300" alt="fstic">
  </a>
</p>

# Fstic: Allele-Frequency-based Genetic distance calculator

[![Documentation](https://img.shields.io/badge/docs-pathogenomics--lab.github.io%2Ffstic-%23d9862b)](https://pathogenomics-lab.github.io/fstic/)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/PathoGenOmics-Lab/fstic/blob/main/LICENSE)
[![fstic](https://img.shields.io/badge/fstic-rust-%23ff8000)](https://github.com/PathoGenOmics-Lab/fstic)
[![Anaconda-Server Badge](https://img.shields.io/conda/dn/bioconda/fstic.svg)](https://anaconda.org/bioconda/fstic)
[![Anaconda-Version Badge](https://anaconda.org/bioconda/fstic/badges/version.svg)](https://anaconda.org/bioconda/fstic)
[![PGO](https://img.shields.io/badge/PathoGenOmics-lab-red?)](https://github.com/PathoGenOmics-Lab)
[![DOI](https://img.shields.io/badge/doi-10.5281%2Fzenodo.16813662-%23ff0077)](https://doi.org/10.5281/zenodo.16813662)

__Paula Ruiz-Rodriguez<sup>1</sup>__
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. I<sup>2</sup>SysBio, University of Valencia-CSIC, FISABIO Joint Research Unit Infection and Public Health, Valencia, Spain </sub>

---

Pairwise genetic distance matrices from allele-frequency data. Fstic reads
single-sample VCFs or frequency tables and computes eight standard estimators
(FST, GST, Jost's *D*, Reynolds, Nei, Cavalli-Sforza chord, Rogers,
Bray-Curtis) in parallel, writing an *N* x *N* matrix.

It is built for data where a sample is not a single genotype but a frequency
spectrum: mixed infections, deep-sequenced isolates, pooled populations.

## 📖 Documentation

**[pathogenomics-lab.github.io/fstic](https://pathogenomics-lab.github.io/fstic/)**

| | |
| --- | --- |
| [**Tutorial**](https://pathogenomics-lab.github.io/fstic/tutorial/) | Four samples from raw VCFs to a tree, finding a mixed infection on the way. Ten minutes, no downloads. |
| [**Installation**](https://pathogenomics-lab.github.io/fstic/getting-started/installation/) | Bioconda, or building from source. |
| [**Input formats**](https://pathogenomics-lab.github.io/fstic/getting-started/inputs/) | What Fstic reads from a VCF or a table, and what it skips. |
| [**Choosing an estimator**](https://pathogenomics-lab.github.io/fstic/guide/metrics/) | All eight formulas, their ranges, and when disagreement between them means something. |
| [**Filtering**](https://pathogenomics-lab.github.io/fstic/guide/filtering/) | Setting thresholds for consensus, within-host and deep amplicon work. |
| [**Command-line reference**](https://pathogenomics-lab.github.io/fstic/guide/cli/) | Every flag, with defaults and exit codes. |
| [**How it works**](https://pathogenomics-lab.github.io/fstic/how-it-works/model/) | The model the estimators rest on, and what stops a run rather than warning. |

## Install

```bash
conda install -c bioconda fstic
```

Or build it, which needs a Rust toolchain:

```bash
git clone https://github.com/PathoGenOmics-Lab/fstic.git
cd fstic && cargo build --release
```

## Sixty seconds

One VCF per sample, named after the sample:

```bash
fstic --vcf TB-*.vcf --output distances.csv
```

```csv
sample,TB-001,TB-002,TB-003
TB-001,0.0000000000,1.0000000000,7.0000000000
TB-002,1.0000000000,0.0000000000,8.0000000000
TB-003,7.0000000000,8.0000000000,0.0000000000
```

Those are cumulative totals, not values in [0, 1]. Adding `--normalize` divides
by the locus count and gives 0.125, 0.875 and 1.0 instead.

Pick an estimator, normalise per locus, and tighten the filters:

```bash
fstic --vcf-list samples.txt \
      --output chord.tsv \
      --formula chord \
      --min-depth 50 \
      --min-af 0.01 \
      --pass-only
```

> [!NOTE]
> `--formula fst` sums the per-locus Nei *G*<sub>ST</sub> across sites, so it
> grows with the locus count and is not bounded by 1. Use `--normalize` for the
> mean, or `--formula gst` for a bounded ratio. The
> [estimator guide](https://pathogenomics-lab.github.io/fstic/guide/metrics/)
> explains the difference.

## Citation

> Ruiz-Rodriguez, P. and Coscolla, M. *Fstic: allele-frequency-based genetic
> distance calculator*. Zenodo.
> [doi:10.5281/zenodo.16813662](https://doi.org/10.5281/zenodo.16813662)

Please cite the original paper for whichever estimator you used as well; they
are listed on the [citation page](https://pathogenomics-lab.github.io/fstic/about/citation/).

## Contributing

Bug reports with a reproducer are the most useful thing you can send. See
[contributing](https://pathogenomics-lab.github.io/fstic/about/contributing/)
for building, testing and adding an estimator.

```bash
cargo test && cargo clippy --all-targets
```


---
<h2 id="contributors" align="center">

✨ [Contributors](https://github.com/PathoGenOmics-Lab/fstic/graphs/contributors)
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



## License

[GPL-3.0](LICENSE)
