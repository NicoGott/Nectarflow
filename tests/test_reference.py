#!/usr/bin/env python3
#
# Tests for the reference library build. Run: python -m unittest discover tests
#
# Three of these are the derivations the coverage figures rest on, so they are
# re-run rather than trusted: the Plant Atlas filter chain, the name
# normalisation rule, and the manifest's own label hygiene. The last one pins
# the rule assign_taxonomy.py has to follow for genus-rank rows: a blank
# species column must never win a species-rank call.
#
import csv, hashlib, os, shutil, subprocess, sys, tempfile, unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "bin"))
import build_reference as BR
import plant_atlas_taxa as PA
import refnames as RN
import resolve_reference as RR

ATLAS = os.path.join(ROOT, "assets", "plant_atlas_2020_gb.tsv")
MANIFEST = os.path.join(ROOT, "assets", "reference_manifest.tsv")
ZIP = os.environ.get("PA2020_ZIP", "")


class Names(unittest.TestCase):
    #
    # The normalisation rule, written down once so a later reader can change it
    # and see exactly what moves
    def test_normalise(self):
        self.assertEqual(RN.normalise("Abelia × grandiflora"), "abelia grandiflora")
        self.assertEqual(RN.normalise("Abelia x grandiflora"), "abelia grandiflora")
        self.assertEqual(RN.normalise("Rosa 'Nevada'"), "rosa")
        self.assertEqual(RN.normalise("Abeliophyllum distichum Roseum Group"), "abeliophyllum distichum roseum")
        self.assertEqual(RN.normalise("  Quercus   robur "), "quercus robur")

    #
    # What may be written at species rank. Open nomenclature never may: these
    # are the strings that put 'Taraxacum sp.' into a results table as a species
    def test_binomial(self):
        self.assertEqual(RN.binomial("Quercus robur L."), "Quercus robur")
        self.assertEqual(RN.binomial("Quercus × bushii"), "Quercus bushii")
        for bad in ("Taraxacum sp.", "Taraxacum agg.", "Syringa x", "Lavandula x", "Hieracium cf.",
                    "Rosa", "Begonia spp.", "Aria x"):
            self.assertEqual(RN.binomial(bad), "", bad)

    def test_genus(self):
        self.assertEqual(RN.genus("Rosa 'Nevada'"), "Rosa")
        self.assertEqual(RN.genus("× Heucherella alba 'Bridget Bloom'"), "Heucherella")


class PlantAtlas(unittest.TestCase):
    #
    # The committed snapshot is re-derived from the Zenodo zip and must come out
    # byte for byte the same. Set PA2020_ZIP to the downloaded file to run it
    @unittest.skipUnless(ZIP and os.path.exists(ZIP), "set PA2020_ZIP to PA2020_NativeStatusData.zip")
    def test_rederive(self):
        blob, _ = PA.load(ZIP)
        rows, counts = PA.derive(blob)
        self.assertEqual(sum(counts.values()), PA.ROWS)
        self.assertEqual(counts["absent"], 27)
        self.assertEqual(len(rows), PA.ROWS - counts["absent"])
        have = RR.read_tsv(ATLAS)
        self.assertEqual(len(have), len(rows))
        self.assertEqual(sorted(r["taxon"] + r["ddbid"] for r in have),
                         sorted(r["taxon"] + r["ddbid"] for r in rows))

    #
    # statusGB, not overallStatus: the dataset covers Ireland and the Channel
    # Islands too, so the filter is load-bearing
    @unittest.skipUnless(os.path.exists(ATLAS), "no committed snapshot")
    def test_statuses(self):
        rows = RR.read_tsv(ATLAS)
        self.assertNotIn("absent", {r["status_gb"] for r in rows})
        self.assertEqual(len(rows), 3470)


@unittest.skipUnless(shutil.which("vsearch"), "vsearch not on PATH")
class Trim(unittest.TestCase):
    #
    # A template buried inside a longer record is found on either strand and the
    # coordinates returned cut exactly that template back out. anchor_all does this
    # with vsearch in one batch, so the test hands it one candidate and one template
    # and checks the slice, which is the property the manifest coordinates rest on.
    def _place(self, seqs, templates):
        with tempfile.TemporaryDirectory() as td:
            return RR.anchor_all("rbcL", seqs, templates, shutil.which("vsearch"), td, 0.9)

    def test_anchor(self):
        import random
        rng = random.Random(1)
        ex = "".join(rng.choice("ACGT") for _ in range(500))
        pad = "".join(rng.choice("ACGT") for _ in range(900))
        cand = pad + ex + pad

        hits = self._place({"cand": cand}, {"tpl": ex})
        self.assertIn("cand", hits)
        strand, s, e = hits["cand"]
        self.assertEqual((strand, cand[s:e]), ("+", ex))

        hits = self._place({"cand": RR.rc(cand)}, {"tpl": ex})
        self.assertIn("cand", hits)
        strand, s, e = hits["cand"]
        self.assertEqual(strand, "-")
        self.assertEqual(RR.rc(cand)[s:e], RR.rc(ex))

    def test_no_anchor(self):
        import random
        rng = random.Random(2)
        hits = self._place({"cand": "".join(rng.choice("ACGT") for _ in range(600))},
                           {"tpl": "".join(rng.choice("ACGT") for _ in range(500))})
        self.assertNotIn("cand", hits)


class Build(unittest.TestCase):
    #
    # The manifest coordinates cut the record, and the md5 is checked against
    # the cut sequence, not the record
    def test_cut_and_check(self):
        raw = "AAAACCCCGGGGTTTT"
        row = {"strand": "+", "start": "5", "end": "12", "accession": "X1.1", "marker": "rbcL", "length": "8"}
        s = BR.cut(row, raw)
        self.assertEqual(s, "CCCCGGGG")
        row["md5"] = hashlib.md5(s.encode()).hexdigest()
        self.assertEqual(BR.check(row, s), s)
        row["md5"] = "0" * 32
        with self.assertRaises(SystemExit):
            BR.check(row, s)

    def test_cut_minus_strand(self):
        row = {"strand": "-", "start": "1", "end": "4"}
        self.assertEqual(BR.cut(row, "AAAACCCCGGGGTTTT"), "AAAA")

    #
    # Identical sequences become one subject, and a label the copies disagree
    # about is dropped rather than picked
    def test_labels(self):
        for r in RR.read_tsv(MANIFEST):
            if r["target_rank"] == "genus":
                self.assertEqual(r["species"], "", r["accession"])
            if r["species"]:
                self.assertEqual(RN.binomial(r["species"]), r["species"], r["accession"])
            self.assertTrue(r["lineage"].startswith("Eukaryota") or "Viridiplantae" in r["lineage"],
                            f"{r['accession']} is not a plant: {r['lineage'][:60]}")

    @unittest.skipUnless(os.path.exists(MANIFEST), "no committed manifest")
    def test_pinned(self):
        for r in RR.read_tsv(MANIFEST):
            self.assertIn(".", r["accession"], "accession must carry its version")
            self.assertEqual(len(r["md5"]), 32)
            self.assertEqual(int(r["end"]) - int(r["start"]) + 1, int(r["length"]))


class AssignRule(unittest.TestCase):
    #
    # A genus-rank reference row has an empty species, and a cluster whose top
    # hits are all such rows must come out at genus rank. The naive test
    # len({species}) == 1 passes on {""} and writes an empty taxon at species
    # rank instead. FAILS until bin/assign_taxonomy.py guards the empty string
    def test_blank_species_cannot_win(self):
        d = tempfile.mkdtemp()
        tax = os.path.join(d, "t.tsv")
        with open(tax, "w") as fh:
            fh.write("accession\tspecies\tgenus\tfamily\torder\tlineage\n")
            for a in ("R1.1", "R2.1"):
                fh.write(f"{a}\t\tRosa\tRosaceae\tRosales\tx\n")
        blast = os.path.join(d, "b.tsv")
        with open(blast, "w") as fh:
            for a in ("R1.1", "R2.1"):
                fh.write(f"q1\t{a}\t99.0" + "\t0" * 8 + "\t900.0\n")
        out = os.path.join(d, "o.tsv")
        subprocess.run([sys.executable, os.path.join(ROOT, "bin", "assign_taxonomy.py"),
                        "--blast", blast, "--taxonomy", tax, "--out", out, "--percentile", "0"], check=True)
        got = list(csv.DictReader(open(out), delimiter="\t"))[0]
        self.assertEqual((got["taxon"], got["rank"]), ("Rosa spp.", "genus"))


if __name__ == "__main__":
    unittest.main()
