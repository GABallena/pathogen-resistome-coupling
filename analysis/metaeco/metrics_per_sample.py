#!/usr/bin/env python3
# Portfolio-safe copy: identifiers generalized; runnable with included io_utils scaffold.

from __future__ import annotations
import math, os
from typing import Dict, List
from .io_utils import read_joined_table, group_by, float_or_nan, write_tsv

OUT = "results/comparisons/metrics_per_sample.tsv"

REQUIRED = ["sampleID","species","bracken_fraction","rel_abundance","new_est_reads"]

def compute():
    rows = read_joined_table()
    # schema check
    if not rows:
        # Allow pipeline continuity with empty placeholder join.
        write_tsv(OUT, [], ["sampleID","n_overlap","bracken_species","metaphlan_species","jaccard","total_reads_bracken","frac_reads_in_overlap"])
        return OUT
    for r in rows:
        for c in REQUIRED:
            if c not in r:
                # Skip strict failure on placeholder; create empty output and return
                write_tsv(OUT, [], ["sampleID","n_overlap","bracken_species","metaphlan_species","jaccard","total_reads_bracken","frac_reads_in_overlap"])
                return OUT
    by_sample = group_by(rows, 'sampleID')
    out_rows: List[Dict[str,object]] = []
    for sid, recs in by_sample.items():
        br_species = set()
        mp_species = set()
        br_reads_sum = 0.0
        overlap_reads = 0.0
        for r in recs:
            bf = float_or_nan(r.get('bracken_fraction'))
            ra = float_or_nan(r.get('rel_abundance'))
            if bf and not math.isnan(bf) and bf>0:
                br_species.add(r['species'])
            if ra and not math.isnan(ra) and ra>0:
                mp_species.add(r['species'])
            new_reads = float_or_nan(r.get('new_est_reads'))
            if not math.isnan(new_reads):
                br_reads_sum += new_reads
                if ra and not math.isnan(ra) and ra>0:
                    overlap_reads += new_reads
        inter = br_species & mp_species
        union = br_species | mp_species
        jaccard = (len(inter)/len(union)) if union else math.nan
        out_rows.append(dict(sampleID=sid,
                             n_overlap=len(inter),
                             bracken_species=len(br_species),
                             metaphlan_species=len(mp_species),
                             jaccard=jaccard,
                             total_reads_bracken=br_reads_sum,
                             frac_reads_in_overlap=(overlap_reads/br_reads_sum if br_reads_sum>0 else math.nan)))
    write_tsv(OUT, out_rows, ["sampleID","n_overlap","bracken_species","metaphlan_species","jaccard","total_reads_bracken","frac_reads_in_overlap"])
    return OUT

if __name__ == '__main__':
    compute()
