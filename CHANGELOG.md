# Changelog

## [1.0.1] - 2026-04-04

### Fixed - second review pass

Numbers change for some inputs. The metrics affected are noted per item.

- **Sites absent from both samples of a pair** produced empty frequency maps rather than "both are reference". `heterozygosity({})` returned 1, so chord added sqrt(2) per such locus (two identical samples came out at 1.41, 2.83, ... depending on cohort size) and GST took 1.0 into its denominator (GST 1.0 dropped to 0.5 purely by adding an unrelated sample to the run). Reynolds inherited it, Nei skipped the locus instead of counting it as identity. *Affects: chord, gst, reynolds, nei*
- **Reference bases read from the wrong contig**: the FASTA lookup fell back to whichever contig sorted first whenever the contig name did not match or the position was past the end, with no warning. Measured Bray-Curtis 0.4 against 0.2 on the same data. Single-contig references still resolve any name; with several contigs an unmatched name is now reported and left un-imputed. `get_base_at` removed. *Affects: all, when `--reference` is used*
- **Allele frequencies summing above 1** (split multi-allelic records with independently estimated frequencies) drove `1 - sum(p^2)` negative. One locus gave FST 1.638, Rogers 1.063, Bray-Curtis 1.300, all bounded by 1. Such sites are now rescaled to a distribution and the count reported. *Affects: all*
- **Reynolds used Nei's GST as theta** instead of the Reynolds, Weir & Cockerham (1983) coancestry coefficient. The two agree only at fixed loci, so published values were not comparable with adegenet, Arlequin or GenAlEx. *Affects: reynolds*
- **Results depended on `--workers`**: rayon's reduction tree varies with thread count and floating point addition is not associative, so a 4-core and a 16-core machine produced different files. Reductions now run over fixed 4096-locus chunks combined in index order. *Affects: all, last digits only*
- **Unparseable VCF FORMAT values failed open**: `AD=".,."` sailed past `--min-alt-reads 50`, since a parse failure became `None` and the filters use `is_none_or`. Now counted and reported.
- **FREQ from AD/DP was never range-checked**, so `AD > DP` yielded a frequency above 1. Range validation moved onto the single path both frequency sources feed into. An unparseable explicit FREQ now falls back to AD/DP and is counted.
- **Inputs vanished from the matrix**: a read error truncated the file at that point (a gzipped VCF produced an empty sample, exit 0), a sample whose variants were all filtered lost its row, and two files with the same basename merged into one. Unreadable files, files with no parseable data lines and duplicate sample names are now errors; a sample with nothing left after filtering keeps its row and gets a warning.
- **Malformed FASTA records were merged**: a nameless record or data before the first header leaked its bases into the next contig, shifting every coordinate after it, and duplicate contig names silently overwrote each other. All three are now errors.
- **Table input did not match its own documentation**: case-insensitive headers and percentage frequencies were promised but unimplemented, and both failed silently (a `Sample,Position,...` header dropped every row, then the run died with an unrelated message). Both now work. Duplicate rows keep the first value, matching the VCF reader, which kept the first while tables kept the last. Deserialise errors keep the file path and message. A table with a `chrom` column can no longer be mixed with one without.
- **The inconsistent-REF warning** counted occurrences rather than positions and named alleles picked by HashMap iteration order, so the message varied between runs.

### Fixed - documentation

- README described FST as the ratio-of-sums, which is what `gst` computes; `fst` sums per-locus GST values and is unbounded without `--normalize`. Both rows corrected.
- README claimed "Indels fully supported" and showed one in its example VCF. Indels, multi-allelic records, `ALT=.` and `ALT=*` are skipped.
- Reynolds theta is now defined in the formula table.
- `--min-alt-rev-reads` and `--pass-only` were missing from the CLI reference.

### Testing

- Weak assertions replaced: three tests checked for the substring `0.000000`, which the hard-coded zero diagonal satisfies on any matrix, one checked that a sample name appeared in the header, and `matrix_symmetric` compared a value against itself.
- `freq_from_ad_dp` used `AD=50,50` in both samples, so it could not detect the ref/alt swap listed as the first critical bug below. It now checks an AD-derived frequency against an explicit FREQ.
- New coverage: table input mode (previously none at all), NA output, sample-name quoting, worker-count independence, cohort-independence of pairwise distances, frequency rescaling, the Reynolds estimator at intermediate frequencies, and the malformed-FASTA cases.

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
