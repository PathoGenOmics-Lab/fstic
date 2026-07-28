# Input formats

Fstic takes either single-sample VCFs or allele-frequency tables. The two are
mutually exclusive in one run, and produce the same numbers on equivalent data.

## VCF

One file per sample. **The sample name is the file stem**, not the name in the
`#CHROM` header line, so `runs/2024/TB-041.vcf` becomes the sample `TB-041`.

```vcf
##fileformat=VCFv4.2
#CHROM POS  ID REF ALT QUAL FILTER INFO FORMAT        sample_A
chr1    100 .  A   T   .    PASS   .    GT:DP:AD:FREQ 0/1:50:24:48.0%
chr1    250 .  C   G   .    PASS   .    GT:DP:AD:FREQ 0/1:45:10:22.2%
```

!!! warning "Two files cannot share a stem"

    `a/s.vcf` and `b/s.vcf` both name the sample `s`. Rather than merge them
    into one row and lose data, Fstic rejects the run. Rename one.

### Fields read

| Field | Required | Used for |
| --- | --- | --- |
| `CHROM`, `POS` | yes | site identity; chromosome is part of the key, so `chr1:100` and `chr2:100` are different loci |
| `REF`, `ALT` | yes | allele identity |
| `FILTER` | no | skipped unless `--pass-only` |
| `FORMAT/FREQ` | preferred | alt allele frequency |
| `FORMAT/DP` | no | `--min-depth` |
| `FORMAT/AD` | no | `--min-alt-reads`, and `FREQ` when it is absent |
| `FORMAT/ADR` | no | `--min-alt-rev-reads` |

Only the first sample column is used. A multi-sample VCF is accepted with a
warning; split it with `bcftools +split` for a real multi-sample analysis.

### Frequency

`FREQ` takes a proportion (`0.48`) or a percentage (`48.0%`). A bare number is
always a proportion.

When `FREQ` is absent, the frequency is computed from `AD/DP`. `AD` is read as
the standard `ref,alt` pair when it carries a comma, and as a bare alt count
otherwise, which is what VarScan writes:

```vcf
chr1  100  .  A  T  .  PASS  .  GT:DP:AD  0/1:100:20,80     # 80 alt of 100 -> 0.8
chr1  100  .  A  T  .  PASS  .  GT:DP:AD  0/1:100:80        # VarScan style -> 0.8
```

A frequency that lands outside `[0, 1]`, which `AD > DP` will do, means the
record cannot be used; those variants are dropped and counted.

### What gets skipped

| Record | Why |
| --- | --- |
| `ALT` containing a comma | multi-allelic; only bi-allelic sites are supported |
| `REF` or `ALT` longer than one base | indels are not supported |
| `ALT` of `.` | monomorphic site |
| `ALT` containing `*` | upstream deletion marker |
| non-`PASS` `FILTER` | only with `--pass-only` |
| fewer than 10 tab-separated columns | no sample column to read, e.g. a sites-only VCF |

Each category is counted and reported at the end of the run.

!!! failure "Uncompressed, tab-delimited only"

    There is no gzip support. A `.vcf.gz` fails with a read error rather than
    producing an empty sample. `bcftools view -Ov` first, or `gunzip`.

## Table

CSV, TSV or TAB, with the delimiter taken from the extension (`.tsv` and `.tab`
are tab-delimited, anything else is comma). Unlike VCF mode, one table can hold
many samples, and `--reference` is required so reference bases can be filled in.

```csv
sample,chrom,position,ref_allele,sequence,frequency,total_dp,alt_dp
sample_A,chr1,100,A,T,0.5,50,24
sample_B,chr1,100,A,T,0.8,60,48
```

| Column | Required | Notes |
| --- | --- | --- |
| `sample` | yes | row and column label in the matrix |
| `position` | yes | 1-based |
| `sequence` | yes | the allele this frequency belongs to |
| `frequency` | yes | proportion (`0.5`) or percentage (`50%`) |
| `chrom` | no | see below |
| `ref_allele` | no | otherwise taken from `--reference` |
| `total_dp` | no | `--min-depth` |
| `alt_dp` | no | `--min-alt-reads` |
| `alt_rv` | no | `--min-alt-rev-reads` |

Column names are case-insensitive and surrounding whitespace is trimmed, so
`Sample,Position,Frequency` works.

### The chrom column

Either every input table carries `chrom` or none does. Mixing the two is
rejected, because the same physical site would be filed under two different
keys and counted twice. With no `chrom` column anywhere, all sites are treated
as belonging to a single unnamed contig.

### Duplicate rows

If the same `(sample, chrom, position, allele)` appears twice, the first value
is kept and the run reports how many duplicates it saw. Contradictory input has
no defensible resolution, so the point is to tell you it is there.

## Reference FASTA

Required for table input, optional for VCF (where `REF` is in the file).
Multi-contig references are supported, and soft-masked (lowercase) sequence is
read as uppercase.

```bash
fstic --table freqs.tsv --reference ref.fa --output d.tsv
```

Contig matching is by name. A **single-contig** reference resolves any name, so
the usual `Chromosome` against `NC_000962.3` mismatch is not a problem. With
several contigs, a name that matches nothing is reported and left un-imputed
rather than silently taking a base from the wrong sequence.

These are rejected outright, since each one shifts coordinates:

- a record with no identifier (`>` on its own)
- sequence data before the first `>` header
- a repeated contig name

## Passing many files

Both modes take a list file instead of a glob:

```bash
fstic --vcf-list paths.txt --output d.csv
fstic --table-list tables.txt --reference ref.fa --output d.csv
```

```text title="paths.txt"
# lines starting with # are ignored
/data/vcfs/TB-001.vcf
/data/vcfs/TB-002.vcf
```
