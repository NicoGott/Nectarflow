#!/usr/bin/env python3
# Swaps rows and columns of a TSV. The ASV tables here are samples-in-rows; mumu wants
# OTUs-in-rows, so the table is turned on its side going in and back again coming out.
import csv
import sys

rows = list(csv.reader(sys.stdin, delimiter="\t"))
w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
w.writerows(zip(*rows))
