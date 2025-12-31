#!/usr/bin/env python3
# Portfolio-safe copy: identifiers generalized; runnable with included io_utils scaffold.

from __future__ import annotations
import math
from typing import Dict, List
from .io_utils import read_joined_table, group_by, float_or_nan, write_tsv, read_tsv

OUT_BRAY = "results/comparisons/beta_bray_curtis_long.tsv"
OUT_AITCHISON = "results/comparisons/beta_aitchison_long.tsv"
CLR_VALUES = "results/comparisons/clr_values.tsv"
FILTER_FLAGS = "results/comparisons/species_filter_flags.tsv"


def bray_curtis(u, v):
    su = sum(u); sv = sum(v)
    if su == 0 and sv == 0: return math.nan
    num = sum(min(a,b) for a,b in zip(u,v))
    den = su + sv
    if den == 0: return math.nan
    return 1 - 2*num/den


def compute():
    joined = read_joined_table()
    # Optional filtering
    keep_species = None
    try:
        flags = read_tsv(FILTER_FLAGS)
        keep_species = {r['species'] for r in flags if r.get('exclude') in ('False','FALSE','0','')}
    except FileNotFoundError:
        pass

    # Build per tool matrices (sample x species)
    by_sample = group_by(joined, 'sampleID')
    species_all_br = []
    species_all_mp = []
    sample_ids = sorted(by_sample.keys())
    for recs in by_sample.values():
        for r in recs:
            sp = r['species']
            if keep_species and sp not in keep_species:
                continue
            if (bf:=float_or_nan(r.get('bracken_fraction'))) and not math.isnan(bf) and bf>0:
                species_all_br.append(sp)
            if (ra:=float_or_nan(r.get('rel_abundance'))) and not math.isnan(ra) and ra>0:
                species_all_mp.append(sp)
    species_br = sorted(set(species_all_br))
    species_mp = sorted(set(species_all_mp))

    def build_matrix(species_list, value_field):
        mat = []
        for sid in sample_ids:
            row_map = {sp:0.0 for sp in species_list}
            for r in by_sample[sid]:
                sp = r['species']
                if keep_species and sp not in species_list: continue
                if sp in row_map:
                    val = float_or_nan(r.get(value_field))
                    if not math.isnan(val):
                        row_map[sp] += val
            mat.append([row_map[sp] for sp in species_list])
        return mat

    br_mat = build_matrix(species_br, 'bracken_fraction') if species_br else []
    mp_mat = build_matrix(species_mp, 'rel_abundance') if species_mp else []

    # Bray-Curtis long
    bray_rows = []
    def pairwise(mat, species_list, tool):
        for i in range(len(sample_ids)):
            for j in range(i+1, len(sample_ids)):
                bc = bray_curtis(mat[i], mat[j])
                bray_rows.append(dict(sample1=sample_ids[i], sample2=sample_ids[j], bray_curtis=bc, tool=tool))
    if len(sample_ids) >= 2:
        if br_mat: pairwise(br_mat, species_br, 'Bracken')
        if mp_mat: pairwise(mp_mat, species_mp, 'MetaPhlAn')
    write_tsv(OUT_BRAY, bray_rows, ["sample1","sample2","bray_curtis","tool"])

    # Aitchison distances from CLR values
    try:
        clr_rows = read_tsv(CLR_VALUES)
    except FileNotFoundError:
        return OUT_BRAY, OUT_AITCHISON
    clr_by_sample = {}
    for r in clr_rows:
        sid = r.get('sampleID'); sp = r.get('species')
        if keep_species and sp not in keep_species: continue
        cb = float_or_nan(r.get('clr_bracken'))
        cr = float_or_nan(r.get('clr_rel'))
        if not math.isnan(cb):
            clr_by_sample.setdefault(('Bracken', sid), {})[sp] = cb
        if not math.isnan(cr):
            clr_by_sample.setdefault(('MetaPhlAn', sid), {})[sp] = cr
    # Align species per tool
    def build_clr_matrix(tool):
        rows = []
        species_all = set()
        for (t,sid), mp in clr_by_sample.items():
            if t==tool:
                species_all.update(mp.keys())
        species_list = sorted(species_all)
        for sid in sample_ids:
            mp = clr_by_sample.get((tool, sid), {})
            row = []
            for sp in species_list:
                val = mp.get(sp, math.nan)
                row.append(val)
            # drop columns with any nan later; for speed just keep
            rows.append(row)
        # drop columns containing any nan
        if not rows: return [], []
        cols_to_keep = [j for j in range(len(rows[0])) if all(not math.isnan(rows[i][j]) for i in range(len(rows)))]
        if not cols_to_keep: return [], []
        matrix = [[rows[i][j] for j in cols_to_keep] for i in range(len(rows))]
        kept_species = [sorted(species_all)[j] for j in cols_to_keep]
        return matrix, kept_species
    import math as _m
    def euclid(u,v): return _m.sqrt(sum((a-b)**2 for a,b in zip(u,v)))
    aitch_rows = []
    for tool in ('Bracken','MetaPhlAn'):
        mat, spp = build_clr_matrix(tool)
        if len(mat) >= 2 and len(spp) >= 2:
            for i in range(len(sample_ids)):
                for j in range(i+1,len(sample_ids)):
                    dist = euclid(mat[i], mat[j])
                    aitch_rows.append(dict(sample1=sample_ids[i], sample2=sample_ids[j], aitchison=dist, tool=tool))
    if aitch_rows:
        write_tsv(OUT_AITCHISON, aitch_rows, ["sample1","sample2","aitchison","tool"])
    else:
        write_tsv(OUT_AITCHISON, [], ["sample1","sample2","aitchison","tool"])
    return OUT_BRAY, OUT_AITCHISON

if __name__ == '__main__':
    compute()
