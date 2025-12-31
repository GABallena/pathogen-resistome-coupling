#!/usr/bin/env python3
# Portfolio-safe copy: identifiers generalized; runnable with included io_utils scaffold.

from __future__ import annotations
import math
from typing import Dict, List
from .io_utils import read_joined_table, group_by, float_or_nan, clr_transform, write_tsv

OUT = "results/comparisons/clr_values.tsv"


def compute():
    rows = read_joined_table()
    by_sample = group_by(rows, 'sampleID')
    combined: Dict[tuple, Dict[str,object]] = {}
    for sid, recs in by_sample.items():
        # collect per tool values
        br_vals = []
        br_species = []
        mp_vals = []
        mp_species = []
        for r in recs:
            bf = float_or_nan(r.get('bracken_fraction'))
            ra = float_or_nan(r.get('rel_abundance'))
            sp = r['species']
            if not math.isnan(bf) and bf>0:
                br_vals.append(bf); br_species.append(sp)
            if not math.isnan(ra) and ra>0:
                mp_vals.append(ra); mp_species.append(sp)
        # transform
        if br_vals:
            br_clr = clr_transform(br_vals)
            for sp, c in zip(br_species, br_clr):
                combined.setdefault((sid, sp), {"sampleID": sid, "species": sp})['clr_bracken'] = c
        if mp_vals:
            mp_clr = clr_transform(mp_vals)
            for sp, c in zip(mp_species, mp_clr):
                combined.setdefault((sid, sp), {"sampleID": sid, "species": sp})['clr_rel'] = c
    out_rows = list(combined.values())
    write_tsv(OUT, out_rows, ["sampleID","species","clr_rel","clr_bracken"])
    return OUT

if __name__ == '__main__':
    compute()
