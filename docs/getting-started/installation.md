# Installation

## Bioconda

The released binary is on Bioconda and is the quickest way in.

=== "conda"

    ```bash
    conda install -c bioconda fstic
    ```

=== "mamba"

    ```bash
    mamba install -c bioconda fstic
    ```

Check it landed:

```bash
fstic --version
```

## From source

Fstic is a single Rust binary with no system dependencies beyond a toolchain.

```bash
git clone https://github.com/PathoGenOmics-Lab/fstic.git
cd fstic
cargo build --release
```

The binary is at `target/release/fstic`. Put it somewhere on your `PATH`:

```bash
install -m 755 target/release/fstic ~/.local/bin/
```

!!! note "Rust version"

    Needs Rust **1.85** or newer, which is checked in CI on every pull request.
    The floor comes from the dependency tree rather than from Fstic's own
    source: `clap_lex` requires edition 2024. On 1.84 the build fails while
    resolving, before it compiles anything.

    `rustup update stable` if `cargo build` complains about `edition2024`.

## Running the test suite

Worth doing once after building from source, particularly on an unusual
platform, since a few of the tests pin exact numeric output.

```bash
cargo test
```

That covers the eight estimators against hand-computed values, the VCF and
table readers, the FASTA reader, output formatting, and a check that results do
not change with the thread count.
