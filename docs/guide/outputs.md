# Output

## The matrix

A full symmetric *N* x *N* matrix with a zero diagonal, a `sample` header
column, and rows in sorted sample order.

```csv
sample,TB-001,TB-002,TB-003
TB-001,0.0000000000,0.3821004184,0.4109589041
TB-002,0.3821004184,0.0000000000,0.2247191011
TB-003,0.4109589041,0.2247191011,0.0000000000
```

Both triangles are written rather than just one, so the file loads directly as
a distance matrix without reshaping.

## Delimiter

Taken from the output extension: `.tsv` and `.tab` give tab-delimited, anything
else gives comma.

```bash
fstic --vcf vcfs/*.vcf -o d.tsv    # tab
fstic --vcf vcfs/*.vcf -o d.csv    # comma
```

## Precision

Ten decimal places. That is more than presentation needs, and deliberately so:
closely related samples in a low-diversity organism produce genuinely small
distances, and at six places an *M. tuberculosis* FST around 1e-8 rounds to
zero, collapsing every near-identical pair into a tie.

## NA cells

Nei's *D* and Reynolds return infinity for a completely fixed difference. That
reaches the file as `NA`, not `inf`:

```csv
sample,fixA,fixB
fixA,0.0000000000,NA
fixB,NA,0.0000000000
```

`NA` is what `read.csv` and `pandas.read_csv` both parse as missing, whereas
`inf` needs special handling in R. If you are getting `NA` where you expected a
number, either the pair genuinely shares no alleles at any locus, or the locus
set is too small for the estimator to say anything.

## Sample names with separators

A name containing the output delimiter, a quote or a newline is quoted per
RFC 4180, with embedded quotes doubled:

```csv
sample,"with,comma",plain
"with,comma",0.0000000000,0.8000000000
plain,0.8000000000,0.0000000000
```

So the column count stays right even with awkward file names.

## stdout and stderr

The matrix goes to the file given by `--output`. Everything else, the filter
echo, warnings, the progress bar and the summary, goes to **stderr**, so
redirecting or piping stdout never mixes diagnostics into data.

To keep a record of a run:

```bash
fstic --vcf-list samples.txt -o d.csv --formula gst 2> run.log
```

## Loading it

=== "R"

    ```r
    m <- as.matrix(read.csv("distances.csv", row.names = 1, check.names = FALSE))
    d <- as.dist(m)

    plot(hclust(d, method = "average"))

    # Neighbour-joining, for chord or Reynolds distances
    library(ape)
    plot(nj(d))
    ```

=== "Python"

    ```python
    import pandas as pd
    from scipy.spatial.distance import squareform
    from scipy.cluster.hierarchy import linkage, dendrogram

    m = pd.read_csv("distances.csv", index_col=0)
    d = squareform(m.values, checks=False)
    dendrogram(linkage(d, method="average"), labels=m.index.tolist())
    ```

=== "scikit-bio"

    ```python
    from skbio import DistanceMatrix
    dm = DistanceMatrix.read("distances.tsv")
    ```

!!! note "NA and downstream tools"

    Most clustering routines will not accept a matrix containing `NA`. If Nei
    or Reynolds gives you some, either drop the affected samples or use a
    bounded estimator such as `gst` or `rogers`, which cannot produce one.
