# Filtering

Five thresholds, applied per variant before anything is computed. The values
Fstic actually used are echoed at the top of every run, and everything dropped
is counted and reported at the end.

## The thresholds

| Flag | Default | Drops a variant when |
| --- | --- | --- |
| `--min-depth` | 30 | `DP` (or `total_dp`) is below it |
| `--min-af` | 0.05 | the alt allele frequency is below it |
| `--min-alt-reads` | 2 | `AD` (or `alt_dp`) is below it |
| `--min-alt-rev-reads` | 2 | `ADR` (or `alt_rv`) is below it |
| `--pass-only` | off | the VCF `FILTER` column is neither `PASS` nor `.` |

```bash
fstic --vcf-list samples.txt --output d.csv \
      --min-depth 50 --min-af 0.01 --min-alt-reads 5 --pass-only
```

## Missing fields do not filter

A threshold can only be applied to a field that is present. If a VCF has no
`ADR`, `--min-alt-rev-reads` cannot do anything, and the variant is kept:

```text
Warning: VCF file 'TB-041.vcf' is missing FORMAT field 'ADR'.
Filtering on this field will be skipped.
```

That warning appears once per file per field. **It is worth reading.** A run
that looks stringent on the command line can be doing far less than you think
if half the fields are absent.

A field that is present but unparseable, such as `AD=".,."`, is a different
case: the filter still cannot be applied, but the run reports how many times it
happened rather than passing over it quietly.

```text
Warning: 1423 FORMAT value(s) for DP/AD/ADR could not be parsed as integers.
The corresponding depth filters could not be applied to those variants.
```

## Frequency has to be usable

A variant with no determinable allele frequency cannot contribute to a
distance, so it is dropped. That covers a missing `FREQ` with no `AD`/`DP` to
fall back on, an unparseable `FREQ` with the same, and any frequency that lands
outside `[0, 1]`, which `AD > DP` produces:

```text
Warning: dropped 12 variant(s) whose allele frequency fell outside [0,1]
(a FREQ out of range, or AD greater than DP).
```

## Setting the frequency threshold

`--min-af` is the one to think hardest about, because it decides what counts as
a real minor allele rather than sequencing noise.

- **Consensus-level analysis**: `--min-af 0.5` and up keeps only the dominant
  allele, so each sample behaves like a single genotype.
- **Within-host diversity**: `0.01` to `0.05`, paired with a real `--min-depth`
  and `--min-alt-reads`. A 1% threshold at 30x depth is meaningless; a single
  read passes it.
- **Deep amplicon data**: below `0.01` only with depth in the thousands and
  strand support enforced through `--min-alt-rev-reads`.

Frequency and depth work together. `--min-alt-reads` is the one that actually
protects against a low-frequency call resting on one or two reads, and
`--min-alt-rev-reads` against a strand artefact.

## What is skipped regardless

Some records are not filtered so much as unsupported, and are dropped whatever
the thresholds say:

- multi-allelic records, meaning `ALT` with a comma
- indels, meaning `REF` or `ALT` longer than one base
- `ALT` of `.` (monomorphic) or containing `*` (upstream deletion)

Each is counted:

```text
Warning: skipped 8814 multi-allelic site(s) (>1 ALT allele).
Only bi-allelic sites are currently supported.
```

If that count is a large fraction of your data, split the records first, which
also gives Fstic a frequency per allele it can use:

```bash
bcftools norm -m-any input.vcf -o split.vcf
```

## Sites summing above 1

Splitting multi-allelic records gives each allele an independently estimated
frequency, and those can add up to more than 1 at a site. Fstic rescales such
sites back to a distribution and says so:

```text
Warning: rescaled 37 site(s) whose allele frequencies summed above 1
(independently estimated frequencies from split multi-allelic records).
```

Without that, `1 - sum(p^2)` goes negative and the estimators leave their
theoretical range. A handful of sites is normal; thousands suggests the caller
is producing frequencies that do not mean what you think.

## Samples with nothing left

A sample can survive filtering with no variants at all. That is legitimate: it
means the sample matches the reference everywhere it was called. It keeps its
row in the matrix and gets a warning:

```text
Warning: TB-093.vcf has no variants left after filtering;
'TB-093' is treated as identical to the reference.
```

If every sample lands there, the run stops with `No polymorphic sites found
after filtering`, which almost always means the thresholds are too strict for
the data or a field name does not match.

## Reading the summary

Every run ends with a block worth checking before the matrix:

```text
--- Summary ---
> Samples: 48
> Loci: 12043
> Pairs: 1128
> Formula: Fst
> Output: distances.csv
> Time: 3.41s
---------------
```

`Loci` is the union of variant sites across all samples after filtering. If it
is far lower than expected, the warnings above the summary will say where the
sites went.
