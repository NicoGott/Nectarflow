#!/usr/bin/env python3
"""Joins the two arms read by read and writes explore/data/crosstab.csv.

Not part of the pipeline. This exists only to build the README figure, and nothing
in replica.nf or reanalysis.nf calls it.

Every raw read gets one of the five fates in figure.py on each side, so the figure's
middle ribbons can be real flows instead of a pairing. SRA spot numbers run densely
from 1, so an accession's raw reads are just 1..N and a fate per read is a byte
array of that length -- the whole join is 80 MB of memory, no read ids held at all.

Arm 1 comes from arm1_read_fate.py. Arm 2 comes from arm2_read_fate.R for everything
from quality filtering onwards, plus the cutadapt outputs still in the work tree to
separate the reads no primer claimed from the reads quality filtering dropped.
"""
import collections, csv, glob, gzip, multiprocessing as mp, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))          # explore/trace
OUT = os.path.join(os.path.dirname(HERE), "data")
# Point REPLICA_RUN at the directory the two pipelines published into: the one
# holding faithful/ and reanalysis/. REPLICA_TRACE is where the per-read fate
# files go, and defaults to a trace/ beside it.
RUN = os.environ.get("REPLICA_RUN") or sys.exit("set REPLICA_RUN to the results directory")
A1, A2 = f"{RUN}/faithful", f"{RUN}/reanalysis"
TRACE = os.environ.get("REPLICA_TRACE", f"{RUN}/trace")
T1, T2 = f"{TRACE}/arm1", f"{TRACE}/arm2"
MARKERS = ("rbcL", "ITS2")

NAMED, QUALITY, PRIMER, SEQUENCE, NAMING = range(5)
FATE = ["named", "quality", "primer", "sequence", "naming"]
# the ranks the assignment step refuses to call a plant; everything else is a name
UNNAMED_RANK = {"unassigned", "order_spanning", "asv_too_short"}


def fasta(path):
    out, h, buf = [], None, []
    for l in open(path):
        if l.startswith(">"):
            if h is not None:
                out.append((h, "".join(buf)))
            h, buf = l[1:].strip(), []
        else:
            buf.append(l.strip())
    if h is not None:
        out.append((h, "".join(buf)))
    return out


def spots(path):
    """-> the SRA spot numbers of a FASTQ, read ids being SRRxxxxxxx.N."""
    with gzip.open(path, "rt") as fh:
        return [int(l.split(".", 1)[1].split()[0])
                for i, l in enumerate(fh) if i % 4 == 0]


# ------------------------------------------------- arm 2, reads cutadapt kept
INDEX = {}


def cutadapt_job(job):
    """One (sample, marker): the spots cutadapt wrote, as a bitset over 1..n_raw."""
    sample, marker, want, n_raw, cands = job
    if want == 0:
        return sample, marker, bytes((n_raw + 7) // 8)
    for p in cands:
        try:
            sp = spots(p)
        except Exception:
            continue
        if len(sp) != want:
            continue
        bits = bytearray((n_raw + 7) // 8)
        for s in sp:
            bits[(s - 1) >> 3] |= 1 << ((s - 1) & 7)
        return sample, marker, bytes(bits)
    return sample, marker, None


def kept_by_cutadapt(sheet, track, threads):
    """-> sample -> bitset of the raw reads that still carried a primer."""
    jobs = []
    for (sample, marker), want in track.items():
        cands = [p for p in sorted(INDEX.get(f"{sample}.{marker}.R1.fastq.gz", ()))
                 if os.path.isfile(p)]
        jobs.append((sample, marker, want, sheet[sample][1], cands))
    out = {s: bytearray((n + 7) // 8) for s, (_a, n) in sheet.items()}
    with mp.Pool(threads) as pool:
        for i, (sample, marker, bits) in enumerate(
                pool.imap_unordered(cutadapt_job, jobs, chunksize=4), 1):
            if bits is None:
                sys.exit(f"no cutadapt output matching the run for {sample} {marker}")
            dst = out[sample]
            for j, b in enumerate(bits):
                dst[j] |= b
            if i % 100 == 0:
                print(f"  cutadapt {i}/{len(jobs)}", flush=True)
    return out


# ----------------------------------------------- arm 2, what a merged read bought
def asv_fate():
    """-> merged sequence -> NAMED or NAMING, for the sequences that survived the
    bimera model. A merged sequence that is not in here was called a chimera."""
    seq2fate = {}
    for m in MARKERS:
        parent = {r["asv"]: r["merged_into"] for r in
                  csv.DictReader(open(f"{A2}/02_asv/{m}.lulu_merges.tsv"), delimiter="\t")}
        rank = {r["asv"]: r["rank"] for r in
                csv.DictReader(open(f"{A2}/03_assignment/{m}.asv_taxonomy.tsv"), delimiter="\t")}
        for h, seq in fasta(f"{A2}/02_asv/{m}.asvs.fasta"):
            a = h.split()[0]
            while a in parent:          # mumu folds an ASV into a retained one
                a = parent[a]
            seq2fate[seq] = NAMING if rank[a] in UNNAMED_RANK else NAMED
    return seq2fate


def main():
    threads = int(os.environ.get("THREADS", 6))
    sheet = {r["sample"]: (r["run"], int(r["read_count"])) for r in
             csv.DictReader(open(f"{RUN}/samplesheet_full.csv"))}
    track = {(r["sample"], r["marker"]): int(r["cutadapt_pairs"]) for m in MARKERS
             for r in csv.DictReader(open(f"{A2}/02_asv/{m}.dada2_track.tsv"), delimiter="\t")
             if r["marker"] == m}

    for root, _d, files in os.walk(f"{RUN}/work"):
        for f in files:
            if f.endswith(".R1.fastq.gz"):
                INDEX.setdefault(f, []).append(os.path.join(root, f))
    print(f"indexed {sum(len(v) for v in INDEX.values()):,} cutadapt outputs", flush=True)

    # one byte per raw read per arm, indexed by spot - 1
    f1 = {s: bytearray([QUALITY]) * n for s, (_a, n) in sheet.items()}
    f2 = {s: bytearray([PRIMER]) * n for s, (_a, n) in sheet.items()}

    for s, bits in kept_by_cutadapt(sheet, track, threads).items():
        arr = f2[s]
        for i in range(len(arr)):
            if bits[i >> 3] >> (i & 7) & 1:
                arr[i] = QUALITY          # cutadapt kept it, so quality filtering lost it
    print("arm 2 primer step done", flush=True)

    CODE = {"named": NAMED, "quality": QUALITY, "primer": PRIMER,
            "sequence": SEQUENCE, "naming": NAMING}
    for p in sorted(glob.glob(f"{T1}/*.tsv.gz")):
        s = os.path.basename(p)[:-7]
        arr = f1[s]
        with gzip.open(p, "rt") as fh:
            for line in fh:
                rid, fate = line.rstrip("\n").split("\t")
                arr[int(rid.split(".", 1)[1]) - 1] = CODE[fate]
    print("arm 1 loaded", flush=True)

    seq2fate = asv_fate()
    print(f"non-chimeric sequences: {len(seq2fate):,}", flush=True)
    for p in sorted(glob.glob(f"{T2}/*.read_fate.tsv.gz")):
        with gzip.open(p, "rt") as fh:
            next(fh)
            for line in fh:
                t = line.rstrip("\n").split("\t")
                arr = f2[t[1]]
                i = int(t[0].split(".", 1)[1]) - 1
                if t[2] != "merged":
                    arr[i] = SEQUENCE     # denoised away, or the pair would not merge
                else:
                    arr[i] = seq2fate.get(t[3], SEQUENCE)   # missing means chimeric
        print(f"  {os.path.basename(p)} loaded", flush=True)

    tab = collections.Counter()
    for s in sheet:
        tab.update(zip(f1[s], f2[s]))
    os.makedirs(OUT, exist_ok=True)
    with open(f"{OUT}/crosstab.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["arm1", "arm2", "reads"])
        for (x, y), n in sorted(tab.items()):
            w.writerow([FATE[x], FATE[y], n])
    print(f"crosstab.csv  {len(tab)} cells, {sum(tab.values()):,} reads", flush=True)


if __name__ == "__main__":
    main()
