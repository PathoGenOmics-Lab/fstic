# Choosing an estimator

Eight estimators, selected with `--formula`. They answer different questions,
and three of them are not bounded by 1, so it is worth knowing which you asked
for before reading the numbers.

Throughout, \(p_i\) and \(q_i\) are the frequencies of allele \(i\) in the two
samples at one locus, \(l\) indexes loci, \(L\) is the number of loci, \(H_S\)
is the mean within-sample heterozygosity and \(H_T\) the heterozygosity of the
pooled pair.

## At a glance

| `--formula` | Range | Scales with *L* | `--normalize` | Good for |
| --- | --- | --- | --- | --- |
| `gst` | [0, 1] | no | no effect | overall differentiation, the safe default if you want a bounded number |
| `fst` | [0, *L*] | yes | divides by *L* | per-locus differentiation, summed |
| `jost_d` | [0, *L*] | yes | divides by *L* | allelic differentiation, insensitive to within-sample diversity |
| `reynolds` | [0, inf) | no | no effect | drift time between recently diverged populations |
| `nei` | [0, inf) | ratio over all loci | no effect | long time scales |
| `chord` | [0, *L*·√2] | yes | divides by *L* | tree building; satisfies the triangle inequality |
| `rogers` | [0, 1] | mean over all loci | no effect | Euclidean distance, already per-locus |
| `bray-curtis` | [0, *L*] | yes | divides by *L* | comparing frequency profiles directly |

## Differentiation

### GST

\[
G_{ST} = \frac{\sum_l (H_{T,l} - H_{S,l})}{\sum_l H_{T,l}}
\]

Nei's (1973) coefficient as a ratio of sums, so informative loci carry more
weight than near-monomorphic ones. Bounded in [0, 1] without any flags, which
makes it the estimator to reach for when you want a number you can compare
across datasets. `--normalize` does not apply.

### FST

\[
F_{ST} = \sum_l \frac{H_{T,l} - H_{S,l}}{H_{T,l}}
\]

The same quantity computed per locus and summed, which is the default. Note the
consequence: this **grows with the number of loci** and is not itself an
\(F_{ST}\) in the bounded sense. Two samples over 40 000 sites can easily give
9150. Pass `--normalize` for the mean per locus, which is bounded.

!!! note "Not Weir & Cockerham"

    This is the Nei formulation. It is not the Weir & Cockerham (1984) \(\theta\)
    estimator, which corrects for sample size and unequal sampling. If a
    reviewer asks for "Weir and Cockerham FST", this is not it.

Summing per-locus ratios weights every locus equally, including ones with a
tiny \(H_T\) where the ratio is noisy. `gst` avoids that by construction. If the
two disagree sharply on your data, that is usually why.

### Jost's *D*

\[
D = \frac{n}{n-1} \cdot \frac{H_T - H_S}{1 - H_S}, \quad n = 2
\]

Measures the fraction of **allelic** diversity partitioned between the two
samples. Unlike \(G_{ST}\) it does not shrink when within-sample diversity is
high, which matters for deeply mixed samples where \(G_{ST}\) can look small
simply because both samples are individually diverse. Summed per locus, so
`--normalize` applies.

## Divergence

### Reynolds

\[
D_R = -\ln(1 - \theta), \quad
\theta = \frac{\sum_l \sum_i (p_i - q_i)^2 / 2}{\sum_l \left(1 - \sum_i p_i q_i\right)}
\]

The Reynolds, Weir & Cockerham (1983) coancestry coefficient, transformed so
that it is linear in drift time for recently diverged populations. This is the
estimator that adegenet, Arlequin and GenAlEx publish under the same name.

Returns infinity, written as `NA`, when \(\theta \ge 1\).

!!! warning "Changed in 1.0.1"

    Earlier versions substituted Nei's \(G_{ST}\) for \(\theta\). The two agree
    only at completely fixed loci, so Reynolds values from before 1.0.1 are not
    comparable with these. For \(p = \{T: 0.5, A: 0.5\}\) against
    \(q = \{T: 1.0\}\) the estimator gives 0.6931 where the substitution gave
    0.4055.

### Nei's standard genetic distance

\[
D = -\ln \frac{J_{xy}}{\sqrt{J_x J_y}}, \quad
J_{xy} = \sum_l \sum_i p_i q_i
\]

Genetic identity turned into a distance. Effective over long time scales, where
it is roughly linear in divergence time under a drift-mutation model.

\(J\) is summed over all loci, including those where both samples are
reference, so shared monomorphic sites pull the distance down. That means Nei's
*D* depends on the size of the locus set, which is the definition working as
intended rather than a quirk. Returns `NA` for a completely fixed difference.

### Cavalli-Sforza chord

\[
D_{CH} = \sum_l \sqrt{2\left(1 - \sum_i \sqrt{p_i q_i}\right)}
\]

A geometric distance on the hypersphere of allele frequencies. It satisfies the
triangle inequality, which is what makes it a reasonable input to
neighbour-joining and other tree builders. Assumes drift without mutation.

Summed over loci, so `--normalize` applies. Per locus it is bounded by √2.

### Rogers

\[
D_R = \sqrt{\frac{\sum_l \sum_i (p_i - q_i)^2}{2L}}
\]

Euclidean distance in frequency space, already divided by the locus count, so
it lands in [0, 1] and `--normalize` does not apply. Being a mean over *L*, it
moves when the locus set grows, the same way Nei's *D* does.

## Profile comparison

### Bray-Curtis

\[
BC = \sum_l \frac{1}{2} \sum_i \lvert p_i - q_i \rvert
\]

A straightforward dissimilarity between frequency profiles, summed over loci.
Per locus it is bounded by 1 and equals the total variation distance.

!!! note "Not a metric"

    Bray-Curtis does not satisfy the triangle inequality, so it should not be
    fed to methods that assume one. Use `chord` for tree building.

## Which to pick

- **A bounded number to compare across datasets**: `gst`, or `fst --normalize`.
- **Building a tree**: `chord`, or `reynolds` for recently diverged samples.
- **Deeply mixed within-host samples**: `jost_d`, since it does not deflate when
  both samples are individually diverse.
- **Deep divergence**: `nei`.
- **Just how different are these frequency vectors**: `bray-curtis` or `rogers`.

Running several is cheap, and disagreement between them is informative rather
than a problem. `gst` far below `jost_d` says within-sample diversity is high;
`fst` far from `gst` says the low-\(H_T\) loci are pulling their weight.
