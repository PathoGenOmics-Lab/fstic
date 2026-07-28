# Command-line reference

```text
fstic [OPTIONS] --output <OUTPUT_FILE>
```

## Input

Exactly one input mode per run. `--vcf` and `--table` cannot be combined.

| Flag | Description |
| --- | --- |
| `--vcf <FILES>...` | One or more single-sample VCF files. Sample name is the file stem. |
| `--vcf-list <FILE>` | A file of VCF paths, one per line. Blank lines and `#` comments ignored. |
| `--table <FILES>...` | One or more frequency tables (`.csv`, `.tsv`, `.tab`). |
| `--table-list <FILE>` | A file of table paths, one per line. |
| `-r`, `--reference <FASTA>` | Reference FASTA. Required for table input, optional for VCF. |

[:octicons-arrow-right-24: Input formats in detail](../getting-started/inputs.md)

## Output

| Flag | Default | Description |
| --- | --- | --- |
| `-o`, `--output <FILE>` | *required* | Destination for the distance matrix. A `.tsv` or `.tab` extension selects tab-delimited output. |

[:octicons-arrow-right-24: Output format](outputs.md)

## Estimator

| Flag | Default | Description |
| --- | --- | --- |
| `-f`, `--formula <NAME>` | `fst` | One of `fst`, `gst`, `nei`, `chord`, `bray-curtis`, `jost_d`, `reynolds`, `rogers`. |
| `--normalize` | off | Divide the cumulative total by the number of loci. Affects `fst`, `chord`, `bray-curtis` and `jost_d`; the others are already ratios or means. |

[:octicons-arrow-right-24: Choosing an estimator](metrics.md)

## Filtering

| Flag | Default | Description |
| --- | --- | --- |
| `--min-depth <INT>` | 30 | Minimum total read depth (`DP`, `total_dp`). |
| `--min-af <FLOAT>` | 0.05 | Minimum alternate allele frequency. |
| `--min-alt-reads <INT>` | 2 | Minimum supporting alt reads (`AD`, `alt_dp`). |
| `--min-alt-rev-reads <INT>` | 2 | Minimum supporting alt reads on the reverse strand (`ADR`, `alt_rv`). |
| `--pass-only` | off | Keep only variants whose VCF `FILTER` is `PASS` or `.`. |

A threshold whose field is absent from the input cannot be applied and the
variant is kept, with a warning per file per field.

[:octicons-arrow-right-24: Filtering in detail](filtering.md)

## Execution

| Flag | Default | Description |
| --- | --- | --- |
| `-w`, `--workers <INT>` | all logical cores | Worker threads. Does not change the result, only the speed. |
| `-h`, `--help` | | Print help. |
| `-V`, `--version` | | Print version. |

## Exit codes

| Code | Meaning |
| --- | --- |
| 0 | Matrix written. |
| 1 | The run stopped; the reason is on stderr. |
| 2 | Argument parsing failed (clap). |

Fstic exits with 1 rather than writing a partial or misleading matrix when:

- an input or reference file cannot be read
- an input VCF yields no parseable data lines, which includes gzipped input
- two inputs share a sample name, since the file stem is the sample name
- a table mixes rows with and without a `chrom` column
- the reference FASTA has a nameless record, data before the first header, or a
  duplicate contig name
- nothing survives filtering, or fewer than two samples remain

## Worked examples

Default FST from a directory of VCFs:

```bash
fstic --vcf vcfs/*.vcf --output distances.csv
```

Bray-Curtis from a list of tables, eight threads:

```bash
fstic --table-list tables.txt \
      --reference ref.fa \
      --output bray_curtis.csv \
      --formula bray-curtis \
      --workers 8
```

Stringent filters, mean per-locus FST, tab-delimited, keeping a log:

```bash
fstic --vcf-list samples.txt \
      --output fst_norm.tsv \
      --formula fst \
      --normalize \
      --min-depth 50 \
      --min-af 0.01 \
      --min-alt-reads 5 \
      --pass-only 2> run.log
```

Several estimators over the same input:

```bash
for f in gst jost_d chord rogers; do
  fstic --vcf-list samples.txt --output "dist_${f}.csv" --formula "$f"
done
```
