#!/usr/bin/env python3
#
# Assigns merged reads to a marker the way the NBGW pipeline does
# (merged_demultiplex_by_primer.py): the forward primer fragment in the first 50
# bases decides, and only when no marker matches there does the
# reverse-complemented reverse fragment in the last 50 get a say.
# Primers are NOT trimmed. Writes one FASTA per marker plus a stats TSV.
#
import argparse, gzip, json, re, sys

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "N": "ACGT"}


#
# Turns a primer pattern into a list of allowed-base sets, one per position
# i.e. 'G[ACT]GACC' -> [{G},{A,C,T},{G},{A},{C},{C}]
def compile_pattern(pattern):
    out = []
    for token in re.findall(r"\[[^\]]+\]|.", pattern):
        out.append(set(token.strip("[]").replace("|", "")) if token.startswith("[") else set(IUPAC.get(token, token)))
    return out


#
# True if the pattern occurs anywhere in the window with at most max_mm
# mismatches. The for/else returns as soon as one position survives the scan
def fuzzy_find(sets, window, max_mm):
    k = len(sets)
    for i in range(len(window) - k + 1):
        mm = 0
        for j in range(k):
            if window[i + j] not in sets[j]:
                mm += 1
                if mm > max_mm:
                    break
        else:
            return True
    return False


#
# Yields (name, sequence) for each record of a FASTQ file, plain or gzipped
def fastq(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline().strip(); fh.readline(); fh.readline()
            yield h[1:].split()[0], s


#
# The original's rule: the first marker whose forward fragment is in the first
# 50 bases wins; if none is, the first whose reverse fragment is in the last 50
def assign(rx, seq, max_mm):
    for m, (start, end, min_len) in rx.items():
        if fuzzy_find(start, seq[:50], max_mm):
            return m
    for m, (start, end, min_len) in rx.items():
        if fuzzy_find(end, seq[-50:], max_mm):
            return m
    return None


#
# Reads the merged reads, gives each one to its marker, and drops those that
# are shorter than the marker's minimum length. Counts the exact matches too,
# so the mismatch tolerance can be judged from the stats file
def main():
    ap = argparse.ArgumentParser(description="Splits merged reads into per-marker FASTA files by their primers")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--merged", required=True, help="merged reads, FASTQ(.gz)")
    ap.add_argument("--markers", required=True,
                    help='JSON: {"rbcL": {"start": regex, "end": regex, "min_len": 450}, ...}')
    ap.add_argument("--mismatches", type=int, default=0,
                    help="mismatches tolerated in each primer fragment (0 = exact, as in the original script)")
    a = ap.parse_args()
    markers = json.loads(a.markers)
    rx = {m: (compile_pattern(v["start"]), compile_pattern(v["end"]), int(v["min_len"])) for m, v in markers.items()}
    out = {m: open(f"{a.sample}.{m}.fasta", "w") for m in rx}
    n = {"merged": 0, "unknown": 0}
    n.update({f"{m}_primer_exact": 0 for m in rx}); n.update({f"{m}_primer": 0 for m in rx})
    n.update({f"{m}_length": 0 for m in rx})
    for name, seq in fastq(a.merged):
        n["merged"] += 1
        # counted but not used to assign: what the original exact-match script would have kept
        exact = assign(rx, seq, 0)
        if exact is not None:
            n[f"{exact}_primer_exact"] += 1
        hit = assign(rx, seq, a.mismatches)
        if hit is None:
            n["unknown"] += 1
            continue
        n[f"{hit}_primer"] += 1
        if len(seq) >= rx[hit][2]:
            n[f"{hit}_length"] += 1
            out[hit].write(f">{name}\n{seq}\n")
    for fh in out.values():
        fh.close()
    with open(f"{a.sample}.split_stats.tsv", "w") as fh:
        fh.write("sample\t" + "\t".join(n) + "\n")
        fh.write(a.sample + "\t" + "\t".join(str(v) for v in n.values()) + "\n")


if __name__ == "__main__":
    sys.exit(main())
