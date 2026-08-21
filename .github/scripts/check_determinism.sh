#!/usr/bin/env bash
# Run every estimator at several thread counts and require byte-identical
# output. A difference here means a reduction is being combined in an order that
# depends on how rayon happened to split the work, which makes the result a
# function of the machine it ran on.
set -euo pipefail

BIN=./target/release/fstic
WORKERS=(1 2 3 5 8)
COMMON=(--vcf det_a.vcf det_b.vcf --min-depth 1 --min-af 0.01
        --min-alt-reads 1 --min-alt-rev-reads 0)

status=0
for formula in fst gst nei chord bray-curtis jost_d reynolds rogers; do
    reference=""
    for w in "${WORKERS[@]}"; do
        out="det_${formula}_w${w}.csv"
        "$BIN" "${COMMON[@]}" --formula "$formula" --workers "$w" --output "$out" 2>/dev/null
        if [[ -z "$reference" ]]; then
            reference="$out"
            continue
        fi
        if ! cmp -s "$reference" "$out"; then
            echo "FAIL ${formula}: --workers ${w} differs from the single-threaded run"
            diff <(cat "$reference") <(cat "$out") || true
            status=1
        fi
    done
    [[ $status -eq 0 ]] && echo "ok   ${formula}"
done

if [[ $status -ne 0 ]]; then
    echo
    echo "Output depends on the thread count. The per-locus reductions in"
    echo "src/calculation/common.rs are meant to run over fixed-size chunks"
    echo "combined in index order; a plain par_iter().sum() reintroduces this."
fi
exit $status
