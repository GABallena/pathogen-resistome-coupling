#!/usr/bin/env python3
# Portfolio-safe copy: identifiers generalized; runnable with included io_utils scaffold.

"""Prepare reduced ordination input matrices (Bray distance & CLR matrix) for later R plotting.
Outputs:
 - comparisons/ordination_bray_bracken.tsv (square matrix)
 - comparisons/ordination_bray_metaphlan.tsv
 - comparisons/ordination_clr_bracken.tsv (samples x features)
 - comparisons/ordination_clr_metaphlan.tsv
"""
from __future__ import annotations
import math
from typing import Dict, List
from .io_utils import read_joined_table, read_tsv, group_by, float_or_nan, write_tsv

FILTER_FLAGS = "results/comparisons/species_filter_flags.tsv"
CLR_VALUES = "results/comparisons/clr_values.tsv"

OUT_BRAY_BR = "results/comparisons/ordination_bray_bracken.tsv"
OUT_BRAY_MP = "results/comparisons/ordination_bray_metaphlan.tsv"
OUT_CLR_BR = "results/comparisons/ordination_clr_bracken.tsv"
OUT_CLR_MP = "results/comparisons/ordination_clr_metaphlan.tsv"


def main():
    joined = read_joined_table()
    keep_species = None
    try:
        flags = read_tsv(FILTER_FLAGS)
        keep_species = {r['species'] for r in flags if r.get('exclude') in ('False','FALSE','0','')}
    except FileNotFoundError:
        pass
    by_sample = group_by(joined, 'sampleID')
    samples = sorted(by_sample.keys())

    def build_matrix(value_col):
        species_all = set()
        for recs in by_sample.values():
            for r in recs:
                if keep_species and r['species'] not in keep_species: continue
                val = float_or_nan(r.get(value_col))
                if not math.isnan(val) and val>0:
                    species_all.add(r['species'])
        species_list = sorted(species_all)
        rows = []
        for sid in samples:
            recs = by_sample[sid]
            mp = {r['species']: float_or_nan(r.get(value_col)) for r in recs}
            row = [mp.get(sp, 0.0) if not math.isnan(mp.get(sp, math.nan)) else 0.0 for sp in species_list]
            rows.append(row)
        return species_list, rows

    def bray_matrix(rows):
        def bray(u,v):
            su=sum(u); sv=sum(v)
            if su==0 and sv==0: return 0.0
            num=sum(min(a,b) for a,b in zip(u,v)); den = su+sv
            return 1 - 2*num/den if den>0 else 0.0
        m = []
        for i in range(len(rows)):
            rline=[]
            for j in range(len(rows)):
                rline.append(bray(rows[i], rows[j]))
            m.append(rline)
        return m

    # Bracken bray
    spp_br, mat_br = build_matrix('bracken_fraction')
    spp_mp, mat_mp = build_matrix('rel_abundance')
    bray_br = bray_matrix(mat_br)
    bray_mp = bray_matrix(mat_mp)

    # write square matrices
    write_tsv(OUT_BRAY_BR, [dict(sample=sid, **{samples[j]: bray_br[i][j] for j in range(len(samples))}) for i,sid in enumerate(samples)], ["sample"]+samples)
    write_tsv(OUT_BRAY_MP, [dict(sample=sid, **{samples[j]: bray_mp[i][j] for j in range(len(samples))}) for i,sid in enumerate(samples)], ["sample"]+samples)

    # CLR matrices (dense wide) for PCA in R
    try:
        clr_rows = read_tsv(CLR_VALUES)
    except FileNotFoundError:
        return
    # partition by tool (non-NA columns)
    clr_by_sample_tool = {('Bracken', sid): {} for sid in samples}
    clr_by_sample_tool.update({('MetaPhlAn', sid): {} for sid in samples})
    for r in clr_rows:
        sid = r.get('sampleID'); sp = r.get('species')
        if keep_species and sp not in keep_species: continue
        cb = float_or_nan(r.get('clr_bracken'))
        cr = float_or_nan(r.get('clr_rel'))
        if not math.isnan(cb): clr_by_sample_tool[('Bracken', sid)][sp] = cb
        if not math.isnan(cr): clr_by_sample_tool[('MetaPhlAn', sid)][sp] = cr
    def write_clr(tool, out_path):
        species_all = set()
        for sid in samples:
            species_all.update(clr_by_sample_tool[(tool,sid)].keys())
        species_list = sorted(species_all)
        rows = []
        for sid in samples:
            mp = clr_by_sample_tool[(tool,sid)]
            row = {"sampleID": sid}
            for sp in species_list:
                row[sp] = mp.get(sp, '')
            rows.append(row)
        write_tsv(out_path, rows, ["sampleID"]+species_list)
    write_clr('Bracken', OUT_CLR_BR)
    write_clr('MetaPhlAn', OUT_CLR_MP)

if __name__ == '__main__':
    main()
