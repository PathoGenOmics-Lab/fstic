# Quick start

## The shortest useful run

One VCF per sample, named after the sample, and an output path:

```bash
fstic --vcf sampleA.vcf sampleB.vcf sampleC.vcf --output distances.csv
```

That computes FST with the default filters. The sample name comes from the file
name, so `sampleA.vcf` becomes the row and column labelled `sampleA`.

The result is a symmetric matrix with a zero diagonal:

<table class="matrix">
<thead><tr><th>sample</th><th>sampleA</th><th>sampleB</th><th>sampleC</th></tr></thead>
<tbody>
<tr><td>sampleA</td><td>0.0000000000</td><td>0.3821004184</td><td>0.4109589041</td></tr>
<tr><td>sampleB</td><td>0.3821004184</td><td>0.0000000000</td><td>0.2247191011</td></tr>
<tr><td>sampleC</td><td>0.4109589041</td><td>0.2247191011</td><td>0.0000000000</td></tr>
</tbody>
</table>

## Many samples

Listing hundreds of files on the command line is awkward, so pass a file of
paths instead. Blank lines and lines starting with `#` are ignored.

```bash
ls vcfs/*.vcf > samples.txt
fstic --vcf-list samples.txt --output distances.csv
```

!!! warning "Paths in a list file"

    Relative paths are resolved against your current directory, not against the
    location of the list file. Use absolute paths, or run from the directory the
    list was written for.

## Choosing an estimator

```bash
fstic --vcf vcfs/*.vcf --output chord.csv --formula chord
```

The eight names are `fst`, `gst`, `nei`, `chord`, `bray-curtis`, `jost_d`,
`reynolds` and `rogers`. [Choosing between them](../guide/metrics.md) covers
what each one measures.

## Per-locus rather than cumulative

FST, chord, Bray-Curtis and Jost's *D* accumulate a value per locus, so they
grow with the number of sites. `--normalize` divides by the locus count to give
the mean instead:

```bash
fstic --vcf vcfs/*.vcf --output fst_mean.tsv --formula fst --normalize
```

A `.tsv` or `.tab` extension switches the output to tab-delimited.

## Tightening the filters

The defaults are deliberately moderate. For a deep-sequenced dataset you will
usually want more:

```bash
fstic \
  --vcf-list samples.txt \
  --output distances.csv \
  --formula fst \
  --min-depth 50 \
  --min-af 0.01 \
  --min-alt-reads 5 \
  --pass-only
```

Everything dropped is counted and reported on stderr at the end of the run, so
check that summary before trusting the matrix. [Filtering](../guide/filtering.md)
goes through each threshold.

## Tables instead of VCFs

If your frequencies already live in a spreadsheet, feed them in directly. Table
input needs `--reference` so the reference base at each site can be recovered.

```bash
fstic --table freqs.csv --reference ref.fa --output distances.csv
```

```csv title="freqs.csv"
sample,chrom,position,ref_allele,sequence,frequency
sampleA,chr1,1043,A,T,0.82
sampleB,chr1,1043,A,T,0.11
```

## Reading the output

The matrix is plain CSV or TSV with a `sample` header column, so it loads
without ceremony:

=== "R"

    ```r
    m <- as.matrix(read.csv("distances.csv", row.names = 1, check.names = FALSE))
    hc <- hclust(as.dist(m), method = "average")
    plot(hc)
    ```

=== "Python"

    ```python
    import pandas as pd
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import squareform

    m = pd.read_csv("distances.csv", index_col=0)
    dendrogram(linkage(squareform(m.values, checks=False), method="average"),
               labels=m.index.tolist())
    ```

Cells that come out non-finite (Nei and Reynolds return infinity for a fixed
difference) are written as `NA`, which both readers understand.
