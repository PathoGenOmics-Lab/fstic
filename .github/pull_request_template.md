## What this changes

<!-- One or two sentences. The commit messages carry the detail. -->

## Why

<!-- The problem, not the patch. For a bug, what went wrong and how it showed. -->

## Does it change any numbers?

<!--
The important box. Several estimators depend on shared code in
src/calculation/common.rs, so a change there moves output for every pair.

If it does: which estimators, and by how much on what input. If it does not,
say so, and say how you know.
-->

- [ ] Output is unchanged for existing inputs
- [ ] Output changes, and the changes are described above and in `CHANGELOG.md`

## Checks

- [ ] `cargo test` passes
- [ ] `cargo clippy --all-targets` adds no warnings
- [ ] `cargo fmt --check` is clean
- [ ] Behaviour changes are covered by a test that fails without the change
- [ ] Documentation under `docs/` moved with the code, if behaviour changed

<!--
On that fourth box: a test that passes both before and after is not covering
the change. Several tests in this repository once asserted that the output
contained "0.000000", which the hard-coded zero diagonal satisfies on any
matrix at all.
-->
