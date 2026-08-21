# Contributing

The full guide lives with the rest of the documentation, at
[pathogenomics-lab.github.io/fstic/about/contributing](https://pathogenomics-lab.github.io/fstic/about/contributing/).
This file is the short version, and the parts specific to working in the
repository itself.

## Reporting a bug

Use the [bug report form](https://github.com/PathoGenOmics-Lab/fstic/issues/new?template=bug_report.yml).
The two things that make a report actionable are the exact command and the
stderr output: Fstic counts everything it drops and prints it before the
summary, so that block usually says what happened.

A wrong number is worth reporting even without a reproducer. Say what you
expected and why: the estimators have variants that differ in ways easy to
mistake for a bug in either direction.

## Before opening a pull request

```bash
cargo test
cargo clippy --all-targets
cargo fmt --check
```

CI runs those, plus a build on the minimum supported Rust version and a check
that the output does not change with the thread count.

## Things worth knowing

**A test that cannot fail is not a test.** Several here once asserted that the
output contained `0.000000`, which the hard-coded zero diagonal satisfies on any
matrix at all. If you add a test with a behaviour change, break the change on
purpose once and check the test notices.

**Identity and fixed-difference cases rarely distinguish two estimators.** Both
send identical samples to 0 and complete differences to 1 or infinity. The case
that tells them apart is an intermediate frequency, which is where the Reynolds
defect hid for a release.

**Use the shared helpers.** `get_all_freqs_at_pos` is where the "an absent
record means reference" rule lives, and `sum_per_locus` is what keeps a
reduction from depending on the thread count. A plain `par_iter().sum()` will
pass the unit tests and fail the determinism job.

**Behaviour changes move with their documentation.** The pages under `docs/`
and the entry in `CHANGELOG.md` are part of the change, not follow-up work. CI
builds the site with `--strict`, so a link to a heading you renamed fails there.

## Mutation testing

Advisory, and run by hand:

```bash
cargo mutants
```

It reports logic no test distinguishes. Configuration is in
`.cargo/mutants.toml`, and the same run is available as a workflow under
Actions.

## Licence

GPL-3.0. Contributions are accepted under the same terms.
