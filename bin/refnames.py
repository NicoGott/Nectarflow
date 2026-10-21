#!/usr/bin/env python3
#
# Plant names: one normalisation rule, one binomial test, and a reader for a
# dated NCBI taxdump. The committed snapshots keep the names exactly as their
# sources wrote them, and every comparison goes through normalise() here, so a
# later reader can change the rule and see precisely what moves.
#
import os, re, sys, unicodedata, zipfile

VIRIDIPLANTAE = 33090
# tokens that mean "this is not a species": a name carrying one is never a binomial
NOT_EPITHET = {"sp", "spp", "cf", "aff", "agg", "x", "group", "hybrid", "var", "subsp", "f", "ssp", "nothosubsp"}


#
# The lookup key for a plant name. Unicode is folded to NFC, the hybrid sign is
# dropped along with a free-standing x, anything inside single quotes is a
# cultivar epithet and goes, a trailing Group goes, whitespace collapses and the
# result is casefolded. 'Abelia x grandiflora' and "Abelia 'Edward Goucher'"
# therefore key as 'abelia grandiflora' and 'abelia'
def normalise(name):
    s = unicodedata.normalize("NFC", str(name)).replace("\xa0", " ")
    s = re.sub(r"[×✕✖]", " ", s)
    s = re.sub(r"'[^']*'|‘[^’]*’", " ", s)
    s = re.sub(r"\b(Group|GROUP)\b", " ", s)
    s = " ".join(t for t in s.split() if t not in ("x", "X"))
    return s.strip().casefold()


#
# The accepted binomial inside a name, or an empty string. Two word tokens, the
# epithet alphabetic and not one of the open-nomenclature tokens above; this is
# what keeps 'Taraxacum sp.', 'Syringa x' and 'Rosa Rugosa Group' from ever
# being written into the library at species rank
def binomial(name):
    t = normalise(name).split()
    if len(t) < 2:
        return ""
    g, e = t[0], t[1].rstrip(".")
    if not g.isalpha() or e in NOT_EPITHET or not re.fullmatch(r"[a-z][a-z-]+", e):
        return ""
    return g.capitalize() + " " + e


#
# The genus inside a name, or an empty string
def genus(name):
    t = normalise(name).split()
    return t[0].capitalize() if t and t[0].isalpha() and len(t[0]) > 2 else ""


class Taxdump:
    #
    # Reads nodes.dmp and names.dmp out of a dated taxdmp zip, keeps only the
    # Viridiplantae subtree, and indexes every name class - not scientific names
    # alone, because 'Mahonia aquifolium' is a synonym and a scientific-name-only
    # lookup silently drops it
    def __init__(self, path):
        self.path = path
        parent, rank = {}, {}
        with zipfile.ZipFile(path) as z:
            with z.open("nodes.dmp") as fh:
                for line in fh:
                    f = line.decode("utf-8", "replace").split("\t|\t")
                    t = int(f[0])
                    parent[t], rank[t] = int(f[1]), f[2]
            keep = _subtree(parent, VIRIDIPLANTAE)
            self.sci, self.index = {}, {}
            with z.open("names.dmp") as fh:
                for line in fh:
                    f = line.decode("utf-8", "replace").split("\t|")
                    t = int(f[0])
                    if t not in keep:
                        continue
                    nm, cls = f[1].strip(), f[3].strip()
                    if cls == "scientific name":
                        self.sci[t] = nm
                    self.index.setdefault(normalise(nm), set()).add(t)
        self.parent, self.rank, self.keep = parent, rank, keep

    #
    # Every Viridiplantae taxid a name maps to, under the normalisation rule
    def lookup(self, name):
        return sorted(self.index.get(normalise(name), ()))

    #
    # The chain of scientific names from the taxon up to Viridiplantae
    def lineage(self, taxid):
        out, t = [], taxid
        while t and t != 1:
            if t in self.sci:
                out.append(self.sci[t])
            if t == VIRIDIPLANTAE:
                break
            t = self.parent.get(t, 0)
        return list(reversed(out))

    #
    # The ancestor at a given rank, by scientific name, or an empty string
    def at_rank(self, taxid, want):
        t = taxid
        while t and t != 1:
            if self.rank.get(t) == want:
                return self.sci.get(t, "")
            t = self.parent.get(t, 0)
        return ""


#
# Every descendant of root, found by inverting the parent map once rather than
# walking to the root for each of two and a half million nodes
def _subtree(parent, root):
    kids = {}
    for t, p in parent.items():
        if t != p:
            kids.setdefault(p, []).append(t)
    seen, stack = {root}, [root]
    while stack:
        for k in kids.get(stack.pop(), ()):
            if k not in seen:
                seen.add(k)
                stack.append(k)
    return seen
