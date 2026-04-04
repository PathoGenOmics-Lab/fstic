# Changelog

## [1.0.1] - 2026-04-04

### Fixed
- **Jost's D formula**: corrected from homozygosity-based to proper heterozygosity-based `D = (n/(n-1)) × (Ht - Hs) / (1 - Hs)` with n/(n-1) correction for 2 populations
- **VCF FREQ parsing**: values without `%` suffix (e.g., `0.5`) are now treated as proportions instead of being divided by 100
- **`--reference` optional for VCF**: reference alleles are already in VCF files; `--reference` is now only required for table inputs
- **Multi-allelic VCF sites**: now warns with count instead of silent discard
- **Empty sample list**: fixed overflow panic when no samples are present
- **Ref allele frequency clamping**: frequencies summing >1.0 due to rounding no longer produce negative reference frequencies

### Changed
- **Deterministic output**: positions sorted before processing; FASTA stored in BTreeMap for deterministic contig order; results reproducible across runs
- **CLI**: migrated from clap 2.34 to clap 4.x with derive macros
- **Error handling**: `main()` returns `ExitCode` via `Result` pattern; no more panics on bad input
- **FASTA reader**: supports multi-contig references stored as `Vec<u8>` (4× memory reduction vs `Vec<char>`)
- **Output format**: auto-detects `.tsv`/`.tab` extension for tab-delimited output; adds `sample` header column
- **Progress**: filter info and progress to stderr; summary with timing at end
- **Thread pool**: gracefully handles double-initialization (no panic in tests)

### Removed
- Redundant `Calculation` enum (uses `Formula` directly)
- Unused type aliases (`Position`, `Allele`, `Frequency`)
- Duplicated `get_all_freqs_at_pos` (was copy-pasted in 7 files, now in `calculation/common.rs`)

### Added
- `calculation/common.rs`: shared `PairSiteFreqs`, `heterozygosity()`, `homozygosity()`, `pooled_heterozygosity()`
- 51 tests (47 unit + 4 integration) covering all metrics, edge cases, symmetry, bounds, FASTA parsing, output format
- Input file existence validation before processing
- `.gitignore`
- Cargo.toml metadata: license, repository, keywords, categories
