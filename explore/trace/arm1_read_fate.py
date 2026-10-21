#!/usr/bin/env python3
"""Per-read fate in arm 1, recovered from files the run already left behind.

Not part of the pipeline. This exists only to build the README figure, and nothing
in replica.nf or reanalysis.nf calls it.

Nothing here re-runs anything. Read ids survive into every arm 1 intermediate, and
vsearch --derep_fulllength is an exact-sequence collapse, so a read's dereplicated
group is recomputable by looking its sequence up in the derep FASTA the run wrote.
From there the cluster map (.uc) and the assignment table finish the chain.

Raw ids are not stored anywhere -- the FASTQs were streamed -- but SRA spot numbers
run densely from 1, so the raw set of an accession is just 1..N for its pair count.
"""
import csv, glob, gzip, os, subprocess, sys

# Point REPLICA_RUN at the directory the two pipelines published into: the one
# holding faithful/ and reanalysis/. REPLICA_TRACE is where the per-read fate
# files go, and defaults to a trace/ beside it.
RUN = os.environ.get("REPLICA_RUN") or sys.exit("set REPLICA_RUN to the results directory")
RES = f"{RUN}/faithful"
OUT = os.environ.get("REPLICA_TRACE", os.path.dirname(RES) + "/trace") + "/arm1"
MARKERS = ("rbcL", "ITS2")


def ids_fastq(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        return [l[1:].split()[0] for i, l in enumerate(fh) if i % 4 == 0]


def fasta(path):
    """-> list of (header, sequence), sequences unwrapped."""
    out, h, buf = [], None, []
    with open(path) as fh:
        for l in fh:
            if l.startswith(">"):
                if h is not None:
                    out.append((h, "".join(buf)))
                h, buf = l[1:].strip(), []
            else:
                buf.append(l.strip())
    if h is not None:
        out.append((h, "".join(buf)))
    return out


INDEX = {}


def build_index():
    """One walk of the work tree. Globbing it per sample is what makes this slow."""
    for root, _dirs, files in os.walk(f"{RUN}/work"):
        for f in files:
            INDEX.setdefault(f, []).append(os.path.join(root, f))
    print(f"indexed {sum(len(v) for v in INDEX.values()):,} work files", flush=True)


def pick(name, want, counter):
    """The work dir holds several attempts; take the one whose size matches the run."""
    for p in sorted(INDEX.get(name, ())):
        if not os.path.isfile(p):      # nextflow stages inputs as links; some targets are gone
            continue
        try:
            if counter(p) == want:
                return p
        except Exception:
            continue
    return None


def load_clusters(marker):
    """derep label -> centroid label, plus the sets of surviving and named centroids."""
    to_centroid = {}
    with open(f"{RES}/05_clusters/{marker}.uc") as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            # the .uc carries the pre-clustering size and the assignment table the summed
            # one, so key on the "sample|index" prefix, which is unique and stable
            if t[0] == "S":
                to_centroid[t[8].split(";")[0]] = t[8].split(";")[0]
            elif t[0] == "H":
                to_centroid[t[8].split(";")[0]] = t[9].split(";")[0]
    # the centroids FASTA is what --minsize 2 left; four clusters in it never reach the
    # assignment table because BLAST returned nothing for them, so take the surviving
    # set from the FASTA and let a missing row mean "no hit", not "dropped"
    present = {h.split(";")[0] for h, _ in fasta(f"{RES}/05_clusters/{marker}.centroids.fasta")}
    named = set()
    with open(f"{RES}/07_assignment/{marker}.assignments.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["rank"] != "discarded":
                named.add(row["cluster"].split(";")[0])
    return to_centroid, present, named


def main():
    os.makedirs(OUT, exist_ok=True)
    track = {r["sample"]: r for r in
             csv.DictReader(open(f"{RES}/read_tracking.tsv"), delimiter="\t")
             if r["sample"] != "TOTAL"}
    only = sys.argv[1:] or sorted(track)
    build_index()
    cl = {m: load_clusters(m) for m in MARKERS}
    print(f"cluster maps loaded: " + ", ".join(f"{m} {len(cl[m][0]):,}" for m in MARKERS), flush=True)

    for s in only:
        t = track[s]
        dst = f"{OUT}/{s}.tsv.gz"
        if os.path.exists(dst):
            continue
        p1 = pick(f"{s}.P1.fastq.gz", int(t["trimmed_pairs"]),
                  lambda p: len(ids_fastq(p)))
        mg = pick(f"{s}.extendedFrags.fastq.gz", int(t["merged"]),
                  lambda p: len(ids_fastq(p)))
        if not p1 or not mg:
            print(f"  {s}: missing trimmed or merged intermediate, skipped", flush=True)
            continue
        fate = {}
        for rid in ids_fastq(p1):
            fate[rid] = "sequence"                      # trimmed but did not merge, unless fixed below
        for rid in ids_fastq(mg):
            fate[rid] = "primer"                        # merged but no marker, unless fixed below
        for m in MARKERS:
            fa = pick(f"{s}.{m}.fasta", int(t[f"{m}_length"]),
                      lambda p: sum(1 for l in open(p) if l[0] == ">"))
            dp = pick(f"{s}.{m}.derep.fasta", -1, lambda p: -1)
            if fa is None or dp is None:
                continue
            seq2lab = {seq: h.split(";")[0] for h, seq in fasta(dp)}
            to_centroid, present, named = cl[m]
            for h, seq in fasta(fa):
                lab = seq2lab.get(seq)
                cen = to_centroid.get(lab) if lab else None
                if cen is None:
                    fate[h.split()[0]] = "sequence"     # its cluster did not survive --minsize 2
                elif cen not in present:
                    fate[h.split()[0]] = "sequence"
                else:
                    fate[h.split()[0]] = "named" if cen in named else "naming"
        with gzip.open(dst, "wt") as out:
            for rid, f in fate.items():
                out.write(f"{rid}\t{f}\n")
        print(f"  {s}: {len(fate):,} reads traced", flush=True)


if __name__ == "__main__":
    main()
