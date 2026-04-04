# Changelog

## [1.0.1] - 2026-04-04

### Fixed — Critical (10 bugs)
- **VCF AD field parsing**: standard comma-separated `AD=ref,alt` was unparseable → `--min-alt-reads` filter never applied
- **Chromosome ignored**: `GenomicPos{chrom, pos}` replaces bare `usize` — chr1:100 and chr2:100 no longer collide
- **FREQ missing → phantom 0.0**: variants without FREQ got freq=0.0 (appeared hom-ref). Now computes AD/DP; skips if neither available
- **"Skipping site" didn't skip**: out-of-range FREQ warning logged but variant still inserted. Now returns None → variant dropped
- **Jost's D formula**: corrected from homozygosity-based to proper heterozygosity-based
- **VCF FREQ parsing**: values without `%` suffix treated as proportions
- **ALT=`.`/`*` processed as alleles**: monomorphic sites and upstream deletions now skipped
- **Indels processed as SNPs**: only single-base REF + single-base ALT accepted (SNPs only)
- **Ref freq overwrite**: `get_all_freqs_at_pos` unconditionally replaced explicit ref freq with imputed value
- **Soft-masked FASTA**: lowercase bases caused allele mismatch (`a` ≠ `A`). Now uppercases on read

### Fixed — Medium (11 bugs)
- **FILTER column ignored**: new `--pass-only` flag skips non-PASS variants
- **CSV malformed rows silent**: now counts and reports skipped rows
- **Table frequencies unvalidated**: values outside [0,1] or NaN/Inf now filtered
- **REF inconsistency across samples**: detects and warns when VCFs used different references
- **Multi-sample VCFs**: warns with count (only first sample column used)
- **Duplicate positions overwrite**: first observation kept via `entry().or_insert()`
- **Sample names with separators**: CSV/TSV output now quotes fields containing separator chars
- **Output precision**: 6 → 10 decimal places (MTB FST can be ~1e-8)
- **FASTA phantom entries**: backfill no longer creates empty entries for every sample×position
- **Infinity as "inf"**: Nei/Reynolds fixed differences now output `NA` for R/Python compatibility
- **`flexible(true)` CSV**: strict column validation; variable-width rows go to error counter
- **`--reference` optional for VCF**: only required for table inputs
- **Empty sample list**: fixed overflow panic
- **Ref allele frequency clamping**: negative frequencies from rounding now clamped to 0

### Changed
- **Deterministic output**: positions sorted before processing; FASTA stored in BTreeMap for deterministic contig order; results reproducible across runs
- **CLI**: migrated from clap 2.34 to clap 4.x with derive macros
- **Error handling**: `main()` returns `ExitCode` via `Result` pattern; no more panics on bad input
- **FASTA reader**: supports multi-contig references stored as `Vec<u8>` (4× memory reduction vs `Vec<char>`)
- **Output format**: auto-detects `.tsv`/`.tab` extension for tab-delimited output; adds `sample` header column
- **Progress**: filter info and progress to stderr; summary with timing at end
- **Thread pool**: gracefully handles double-initialization (no panic in tests)
- **VcfVariant.alt_freq**: now `f64` (not `Option<f64>`) — variants without determinable frequency skipped at parse time

### Removed
- Redundant `Calculation` enum (uses `Formula` directly)
- Unused type aliases (`Position`, `Allele`, `Frequency`)
- Duplicated `get_all_freqs_at_pos` (was copy-pasted in 7 files, now in `calculation/common.rs`)

### Added
- `GenomicPos{chrom, pos}` type for chromosome-aware position tracking
- `--pass-only` flag for VCF FILTER column support
- FREQ computation from AD/DP when FREQ field is absent
- `calculation/common.rs`: shared `PairSiteFreqs`, `heterozygosity()`, `homozygosity()`, `pooled_heterozygosity()`
- 57 tests (49 unit + 8 integration) covering all metrics, edge cases, symmetry, bounds, FASTA parsing, output format, indels, pass-only, freq-from-AD/DP
- Input file existence validation before processing
- `.gitignore`
- Cargo.toml metadata: license, repository, keywords, categories
