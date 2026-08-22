#!/usr/bin/env python3
"""Write two VCFs large enough that the summation order shows in the output.

Floating point addition is not associative, so a reduction whose tree depends on
the thread count differs in the last digits. Making that visible needs care in
two directions, both of which were established by breaking the reduction on
purpose and checking that this input catches it:

  * Size. On a handful of loci the difference sits below the printed precision
    and nothing is observable. The effect this guards against was found at a few
    hundred thousand loci, so that is the scale used here.

  * Irregularity. An earlier version of this generator walked the frequencies
    with a fixed arithmetic stride, and the deliberately broken build passed:
    values that regular produce partial sums that agree whatever order they are
    added in. Drawing them at random is what exposes it.

Seeded, so "random" here means irregular rather than different every run. A
failure reproduces exactly, which is the whole point of a regression check.
"""

import pathlib
import random

POSITIONS = 400_000
DENSITY = 0.85
SEED = 20260728

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=4411532>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
)


def write(name: str, rng: random.Random) -> None:
    rows = [HEADER]
    for pos in range(1, POSITIONS + 1):
        if rng.random() >= DENSITY:
            continue
        freq = rng.uniform(2.0, 98.0)
        rows.append(
            f"chr1\t{pos}\t.\tA\tT\t.\tPASS\t.\tGT:DP:FREQ\t0/1:100:{freq:.6f}%\n"
        )
    pathlib.Path(name).write_text("".join(rows))
    print(f"wrote {name} ({len(rows) - 1} variants)")


rng = random.Random(SEED)
write("det_a.vcf", rng)
write("det_b.vcf", rng)
