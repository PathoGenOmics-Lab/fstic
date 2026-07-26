# Reproducibility and failure modes

## The same input gives the same bytes

Running Fstic twice on the same files produces byte-identical output, including
on a machine with a different core count.

Three things have to hold for that, and each is arranged deliberately:

**Loci and samples are sorted** before the pair loop, so processing order does
not depend on the order files were passed or on hash iteration order.

**Reductions run over fixed-size chunks.** Floating point addition is not
associative, so the order in which per-locus terms are summed shows up in the
last digits. rayon splits an indexed iterator according to the current thread
count and work stealing, which would make the result a function of `--workers`.
Fstic instead sums chunks of 4096 loci sequentially and combines them in index
order, so the summation tree depends on the number of loci and nothing else.

**Duplicate records resolve by position, not by timing.** Where the same
`(sample, locus, allele)` appears twice, the first occurrence wins in both
readers, and the readers preserve input order.

!!! warning "Changed in 1.0.1"

    Before 1.0.1 the reduction was thread-dependent. On 320 000 loci, `-w 3`
    gave FST 91356.2127736697 where every other thread count gave ...699. A
    4-core and a 16-core machine produced different files for the same command.

`--workers` therefore changes only how long the run takes.

## What stops a run

Fstic exits with status 1 rather than writing a matrix it cannot stand behind:

| Condition | Why it is fatal |
| --- | --- |
| An input or reference file cannot be read | nothing to compute from |
| A VCF yields no parseable data lines | includes gzipped input, which fails on line one |
| Two inputs share a sample name | the file stem is the sample name, so they would merge into one row |
| A table mixes rows with and without `chrom` | the same site would be filed under two keys and counted twice |
| A nameless FASTA record, or data before the first header | its bases would leak into the next contig and shift every coordinate |
| A duplicate contig name in the FASTA | which sequence wins would depend on file order |
| Nothing survives filtering, or fewer than two samples | no matrix to write |

The common thread is that each of these has a plausible-looking wrong answer
available, and producing it quietly is worse than stopping.

## What only warns

Conditions that are legitimate, or that lose data in a bounded and reportable
way, warn and continue. All of them are counted and printed before the summary:

- skipped multi-allelic sites, indels, `ALT=.` and `ALT=*`
- `FORMAT` values that could not be parsed, so a filter could not be applied
- variants dropped for a frequency outside `[0, 1]`
- sites rescaled because their frequencies summed above 1
- contig names absent from the reference, and positions past a contig's end
- a sample with no variants left after filtering, which keeps its row
- inconsistent `REF` alleles across samples at the same position
- multi-sample VCFs, where only the first sample column is read
- data lines with fewer than ten columns, e.g. a sites-only VCF

!!! tip "Keep the log"

    Warnings go to stderr, so `2> run.log` captures them without touching the
    matrix. The counts are the quickest way to tell a clean run from one that
    quietly threw away half the data.

## Inconsistent references

If two VCFs were called against different references, the same position can
carry different `REF` alleles. Fstic detects it and names the first position
that disagrees:

```text
Warning: 412 position(s) have inconsistent REF alleles across samples,
first at chr1:1043 (A vs G). Check that all inputs were called against
the same reference.
```

This is a warning rather than an error because a handful of positions can
disagree for benign reasons. Hundreds means the inputs are not comparable and
the matrix is not meaningful.

## Numerical bounds

Each estimator is bounded, and the tests assert it. Two of those bounds depend
on frequency vectors being distributions, which is why the rescaling step
exists:

| Estimator | Bound |
| --- | --- |
| `gst`, `rogers` | [0, 1] |
| `fst`, `jost_d`, `bray-curtis` | [0, *L*], or [0, 1] with `--normalize` |
| `chord` | [0, *L*·√2] |
| `nei`, `reynolds` | [0, inf), written as `NA` when infinite |

A value outside these is a bug. Before 1.0.1 one locus with frequencies summing
to 1.8 was enough to give FST 1.638 and Rogers 1.063.

## Test coverage

`cargo test` runs 63 unit tests and 20 integration tests covering:

- each estimator against hand-computed values, plus symmetry and bounds
- identity of indiscernibles, and that distances do not shift when an unrelated
  sample joins the run
- Reynolds at intermediate frequencies, where the coancestry estimator and
  Nei's GST diverge
- VCF parsing: `AD` in both layouts, frequency from `AD/DP`, percentages,
  indels, `ALT=*`, `--pass-only`, unparseable fields
- table parsing, and agreement with VCF mode on equivalent data
- FASTA parsing, including the three malformed shapes that are rejected
- output formatting: delimiters, `NA`, quoting, precision
- byte-identical output at 1, 2, 3 and 8 workers

The integration tests drive the real binary through `CARGO_BIN_EXE_fstic`, so
they exercise argument parsing and exit codes as well.
