#!/usr/bin/env python3
#
# Cuts the PCR primers off the query sequences before BLAST, using the same
# primer fragments and mismatch tolerance as split_by_primer.py. A sequence
# where either fragment cannot be found goes out unchanged, and is counted.
# Why we trim at all: docs/replica.md#primer-trimming-of-blast-queries
#
import argparse, json, sys
sys.path.insert(0, __import__("os").path.dirname(__file__))
from split_by_primer import compile_pattern


#
# Start index of the best (fewest mismatches) occurrence within max_mm, or None.
# With last=True the window is scanned backwards, so the reverse primer nearest
# the end of the sequence wins
def find(sets, window, max_mm, last=False):
    k, best = len(sets), None
    rng = range(len(window) - k, -1, -1) if last else range(len(window) - k + 1)
    for i in rng:
        mm = sum(1 for j in range(k) if window[i + j] not in sets[j])
        if mm <= max_mm and (best is None or mm < best[1]):
            best = (i, mm)
            if mm == 0:
                break
    return None if best is None else best[0]


#
# Streams the FASTA, trims each sequence and writes it to stdout, then reports
# how many were left untrimmed
def main():
    ap = argparse.ArgumentParser(description="Trims the PCR primers off the BLAST query sequences",
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--marker", required=True)
    ap.add_argument("--markers", required=True, help="same JSON as split_by_primer.py")
    ap.add_argument("--mismatches", type=int, default=2)
    a = ap.parse_args()
    m = json.loads(a.markers)[a.marker]
    start, end = compile_pattern(m["start"]), compile_pattern(m["end"])
    n = untrimmed = 0
    name, seq = None, []

    #
    # Writes the sequence just read, cut between the end of the forward primer
    # and the start of the reverse one. Both are searched in a 50 base window,
    # the same window split_by_primer.py used to recognise the marker
    def emit():
        nonlocal n, untrimmed
        s = "".join(seq); n += 1
        i = find(start, s[:50], a.mismatches)
        j = find(end, s[-50:], a.mismatches, last=True)
        if i is None or j is None:
            untrimmed += 1; t = s
        else:
            t = s[i + len(start): len(s) - 50 + j]
        sys.stdout.write(f">{name}\n{t}\n")

    for line in open(a.fasta):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                emit()
            name, seq = line[1:], []
        else:
            seq.append(line)
    # the last record has no following header to flush it
    if name is not None:
        emit()
    print(f"{a.marker}: {n} queries, {untrimmed} left untrimmed (primer fragment not found)", file=sys.stderr)


if __name__ == "__main__":
    main()
