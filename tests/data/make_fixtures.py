#!/usr/bin/env python3
"""Write the polyscan test fixtures into this directory.

Deterministic (fixed seed, gzip mtime 0), so re-running reproduces the
committed files byte for byte. Stdlib only.

    python3 tests/data/make_fixtures.py

Background sequence is "GC" repeats, which hold no A or T, so every run below is
isolated unless a fixture deliberately places runs near each other. Comments give
the expected tracts at the default W10/P80 (threshold 8 of 10).
"""

import gzip
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent


def bg(n):
    """n bases of A/T-free background."""
    return ("GC" * (n // 2 + 1))[:n]


def write_fasta(name, records, width=60):
    with open(HERE / name, "w") as fh:
        fh.write(fasta_text(records, width))


def fasta_text(records, width=60):
    out = []
    for rid, seq in records:
        out.append(f">{rid}\n")
        for i in range(0, len(seq), width):
            out.append(seq[i:i + width] + "\n")
    return "".join(out)


def write_gz(name, text):
    with open(HERE / name, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
            gz.write(text.encode())


def main():
    # Contigs shorter than w are skipped; an empty FASTA yields no output.
    write_fasta("short_contig.fa", [("short", "AAAAA"), ("also_short", "ACGTACGTA")])
    (HERE / "empty.fa").write_text("")

    # A run touching position 0 and one running to the last base.
    write_fasta("edges.fa", [("edges", "A" * 10 + bg(30) + "A" * 12)])

    # An A8 block gives passing starts s-2..s, i.e. one tract [s-2, s+10).
    # bookended: second block 4 bases later -> tracts touch (end == start), merge.
    # gap1:      second block 5 bases later -> tracts 1 base apart, do not merge.
    # overlap:   second block 3 bases later -> starts 21..28 fail, but tracts
    #            [18,30) and [29,41) overlap, so they must merge.
    write_fasta("bookended.fa", [("bookended", bg(20) + "A" * 8 + bg(4) + "A" * 8 + bg(20))])
    write_fasta("gap1.fa", [("gap1", bg(20) + "A" * 8 + bg(5) + "A" * 8 + bg(20))])
    write_fasta("overlap_gap.fa", [("overlap_gap", bg(20) + "A" * 8 + bg(3) + "A" * 8 + bg(20))])

    # A20 then T10: windows straddling the boundary pass on '+' (8A2T) and '-'
    # (2A8T), so a '-' tract starts before the longer '+' tract ends.
    # AT repeat: at p=50 every window passes on both strands at the same coords.
    write_fasta("both_strands.fa", [
        ("a_then_t", bg(20) + "A" * 20 + "T" * 10 + bg(20)),
        ("at_repeat", bg(20) + "AT" * 10 + bg(20)),
    ])

    # Soft-masked run, an N run, and A runs interrupted by IUPAC codes (which
    # count towards no base).
    write_fasta("softmask_n_iupac.fa", [
        ("softmask", bg(20) + "a" * 10 + bg(20) + "t" * 4 + "T" * 6 + bg(20)),
        ("n_run", bg(20) + "N" * 15 + bg(20)),
        ("iupac", bg(20) + "AAAARAAAAYAA" + bg(20) + "AAAAAAAARY" + bg(20)),
    ])

    # Several contigs in non-alphabetical order, plus a gzipped copy.
    multi = [
        ("chrZ", bg(15) + "A" * 11 + bg(15)),
        ("chr2", bg(15) + "T" * 9 + "G" + "T" * 3 + bg(15)),
        ("chr10", bg(5)),  # shorter than w
        ("chrA", bg(15) + "A" * 10 + bg(10) + "T" * 10 + bg(15)),
    ]
    write_fasta("multi_contig.fa", multi)
    write_gz("multi_contig.fa.gz", fasta_text(multi))

    # Float-rounding cases (see v0.2 plan, Bug fixes):
    # thresh_p56_w25: 14 A in a 25-window. Exact threshold ceil(56*25/100) = 14,
    #                 but v0.1.0 computes 15, so it reports nothing at -w 25 -p 56.
    # score_w20:      11 A in a 20-window. Exact score ceil(100*11/20) = 55,
    #                 but v0.1.0 prints 56 at -w 20 -p 55.
    write_fasta("rounding.fa", [
        ("thresh_p56_w25", bg(30) + "A" * 14 + bg(30)),
        ("score_w20", bg(30) + "A" * 11 + bg(30)),
    ])

    # Seeded random 1 Mb contig with planted A and T runs of varied length and
    # purity. Not used for windows goldens (output would be large); for unit
    # tests and the Phase 2 bedtools-equivalence goldens.
    rng = random.Random(20261018)
    seq = [rng.choice("ACGT") for _ in range(1_000_000)]
    for _ in range(2000):
        base = rng.choice("AT")
        length = rng.randint(6, 40)
        purity = rng.choice([1.0, 0.95, 0.9, 0.8, 0.7])
        start = rng.randrange(0, len(seq) - length)
        for i in range(start, start + length):
            seq[i] = base if rng.random() < purity else rng.choice("ACGT")
    write_gz("random_1mb.fa.gz", fasta_text([("random_1mb", "".join(seq))]))


if __name__ == "__main__":
    main()
