# Contributing

Bug reports, especially ones with a reproducer, are the most useful thing you
can send.

## Reporting a bug

Open an [issue](https://github.com/PathoGenOmics-Lab/fstic/issues) with:

- the exact command
- `fstic --version`
- the stderr output, which carries the filter echo, the warning counts and the
  summary
- a minimal input that shows the problem, ideally a couple of small VCFs

A wrong number is worth reporting even without a reproducer. Say what you
expected and why, since the estimators have variants that differ in ways that
are easy to mistake for a bug and vice versa.

## Building and testing

```bash
git clone https://github.com/PathoGenOmics-Lab/fstic.git
cd fstic
cargo build
cargo test
cargo clippy --all-targets
```

Please keep `cargo test` green and add no new clippy warnings.

## Adding an estimator

Each lives in its own file under `src/calculation/` and exposes one function
taking the two samples' `PositionalData`, the sorted locus list, and whatever
flags it honours.

Use the shared helpers rather than rolling your own loop:

- `get_all_freqs_at_pos` builds both frequency vectors at a locus, with the
  reference allele imputed. It is where the "absent means reference" rule
  lives, so calling it is what keeps a new estimator consistent with the rest.
- `sum_per_locus` and `sum_per_locus_n` accumulate per-locus terms in a
  thread-count-independent order. A plain `par_iter().sum()` will make the
  output depend on `--workers`.
- `heterozygosity` and `pooled_heterozygosity` for the GST family.

Then register it in the `Formula` enum in `src/cli.rs` and the match in
`src/main.rs`.

## Tests a new estimator needs

Look at `src/calculation/tests.rs` for the shape. At minimum:

- identical samples give exactly 0
- symmetry: `d(a, b) == d(b, a)`
- a hand-computed value at intermediate frequencies, which is the one that
  catches a wrong formula. Identity and fixed-difference cases often agree
  between two different estimators, so they cannot tell them apart on their own.
- the theoretical bound
- if the estimator is cumulative, that a locus absent from both samples
  contributes nothing (see `absent_site_is_neutral_for_cumulative_metrics`)

Assert on parsed numbers, not on substrings of the output file. The diagonal is
hard-coded to zero, so `content.contains("0.000000")` passes on any matrix at
all; the `off_diagonal` helper in `tests/integration.rs` exists for this.

## Documentation

The site is MkDocs Material.

```bash
pip install -r requirements-docs.txt
mkdocs serve
```

CI builds with `--strict`, so a broken internal link or a page missing from the
nav fails the build rather than shipping. If you change behaviour, the page
describing it and `CHANGELOG.md` should move with the code.

## Licence

GPL-3.0. Contributions are accepted under the same terms.
