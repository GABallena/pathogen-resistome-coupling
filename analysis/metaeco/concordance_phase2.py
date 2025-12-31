#!/usr/bin/env python3
# Portfolio-safe copy: identifiers generalized; runnable with included io_utils scaffold.

from __future__ import annotations
import math
from typing import List, Dict
from .io_utils import read_joined_table, read_tsv, write_tsv, float_or_nan, group_by, spearman

OUT_SAMPLE = "results/comparisons/sample_concordance_phase2.tsv"
OUT_SPECIES_PREV = "results/comparisons/species_prev_stratified_corr.tsv"
OUT_SPECIES_BA = "results/comparisons/species_bland_altman_clr.tsv"

CLR_VALUES = "results/comparisons/clr_values.tsv"
FILTER_FLAGS = "results/comparisons/species_filter_flags.tsv"

# reuse Lin's CCC + Bland Altman from baseline (could refactor; replicate for clarity)

def lins_ccc(x: List[float], y: List[float]) -> float:
    if len(x) != len(y) or len(x) < 3: return math.nan
    pair = [(a,b) for a,b in zip(x,y) if not(math.isnan(a) or math.isnan(b))]
    if len(pair) < 3: return math.nan
    xs = [p[0] for p in pair]; ys = [p[1] for p in pair]
    mx = sum(xs)/len(xs); my = sum(ys)/len(ys)
    vx = sum((a-mx)**2 for a in xs)/(len(xs)-1)
    vy = sum((b-my)**2 for b in ys)/(len(ys)-1)
    sxy = sum((a-mx)*(b-my) for a,b in pair)/(len(xs)-1)
    denom = vx + vy + (mx - my)**2
    if denom == 0: return math.nan
    return (2*sxy)/denom


def bland_altman(br, mp):
    pair = [(a,b) for a,b in zip(br, mp) if not(math.isnan(a) or math.isnan(b))]
    if len(pair) < 3: return (math.nan, math.nan, math.nan, math.nan)
    diffs = [a-b for a,b in pair]
    bias = sum(diffs)/len(diffs)
    var = sum((d-bias)**2 for d in diffs)/(len(diffs)-1) if len(diffs)>1 else math.nan
    import math as _m
    sd = _m.sqrt(var) if not _m.isnan(var) else math.nan
    loa_lower = bias - 1.96*sd if sd and not _m.isnan(sd) else math.nan
    loa_upper = bias + 1.96*sd if sd and not _m.isnan(sd) else math.nan
    return (bias, loa_lower, loa_upper, sum(abs(d) for d in diffs)/len(diffs))


def compute():
    # load filter flags
    flags = {r['species']: r for r in read_tsv(FILTER_FLAGS)} if FILTER_FLAGS and len(FILTER_FLAGS)>0 else {}
    keep = {sp for sp, fr in flags.items() if fr.get('exclude','FALSE') in ('False','FALSE','0','')}
    # load joined & clr
    joined = read_joined_table()
    clr_rows = read_tsv(CLR_VALUES)
    clr_idx = {(r['sampleID'], r['species']): (float_or_nan(r.get('clr_bracken')), float_or_nan(r.get('clr_rel'))) for r in clr_rows}
    by_sample = group_by(joined, 'sampleID')
    sample_rows = []
    for sid, recs in by_sample.items():
        frac_pairs = []
        clr_pairs = []
        for r in recs:
            if keep and r['species'] not in keep: continue
            bf = float_or_nan(r.get('bracken_fraction'))
            ra = float_or_nan(r.get('rel_abundance'))
            crb, crr = clr_idx.get((sid, r['species']), (math.nan, math.nan))
            frac_pairs.append((bf, ra))
            clr_pairs.append((crb, crr))
        bf_list = [a for a,_ in frac_pairs]; ra_list=[b for _,b in frac_pairs]
        cb_list = [a for a,_ in clr_pairs]; cr_list=[b for _,b in clr_pairs]
        ccc_frac = lins_ccc(bf_list, ra_list)
        ccc_clr = lins_ccc(cb_list, cr_list)
        bias, loa_l, loa_u, madiff = bland_altman(cb_list, cr_list)
        sample_rows.append(dict(sampleID=sid, lins_ccc_fraction=ccc_frac, lins_ccc_clr=ccc_clr,
                                bland_altman_bias_clr=bias, bland_altman_loa_lower=loa_l, bland_altman_loa_upper=loa_u,
                                mean_abs_diff_clr=madiff))
    write_tsv(OUT_SAMPLE, sample_rows, ["sampleID","lins_ccc_fraction","lins_ccc_clr","bland_altman_bias_clr","bland_altman_loa_lower","bland_altman_loa_upper","mean_abs_diff_clr"])

    # prevalence stratified correlations (filtered species)
    by_species = group_by([r for r in joined if (not keep) or r['species'] in keep], 'species')
    by_sample_keys = sorted(by_sample.keys())
    sid_index = {sid:i for i,sid in enumerate(by_sample_keys)}
    species_rows = []
    for sp, recs in by_species.items():
        br_vec=[0.0]*len(by_sample_keys); mp_vec=[0.0]*len(by_sample_keys)
        for r in recs:
            idx = sid_index[r['sampleID']]
            bf = float_or_nan(r.get('bracken_fraction')); ra = float_or_nan(r.get('rel_abundance'))
            if not math.isnan(bf): br_vec[idx]+=bf
            if not math.isnan(ra): mp_vec[idx]+=ra
        nz = [i for i,(a,b) in enumerate(zip(br_vec, mp_vec)) if a>0 or b>0]
        if len(nz) >= 3:
            rho = spearman([br_vec[i] for i in nz],[mp_vec[i] for i in nz])
        else:
            rho = math.nan
        prev_br = sum(1 for a in br_vec if a>0); prev_mp=sum(1 for a in mp_vec if a>0)
        prev_any = max(prev_br, prev_mp)
        bins = [0,1,3,5,10,20,10**9]; labels=["0","1","2-3","4-5","6-10","11-20","21+"]
        for b1,b2,label in zip(bins[:-1], bins[1:], labels):
            if b1 <= prev_any < b2: prev_bin=label; break
        species_rows.append(dict(species=sp, spearman=rho, prev_br=prev_br, prev_mp=prev_mp, prev_bin=prev_bin))
    write_tsv(OUT_SPECIES_PREV, species_rows, ["species","spearman","prev_br","prev_mp","prev_bin"])

    # Bland-Altman CLR per species
    ba_rows=[]
    clr_by_species = group_by([r for r in clr_rows if (not keep) or r['species'] in keep], 'species')
    for sp, recs in clr_by_species.items():
        br_list=[]; mp_list=[]
        for r in recs:
            cb = float_or_nan(r.get('clr_bracken')); cr=float_or_nan(r.get('clr_rel'))
            if not math.isnan(cb) and not math.isnan(cr):
                br_list.append(cb); mp_list.append(cr)
        if not br_list: continue
        bias, loa_l, loa_u, madiff = bland_altman(br_list, mp_list)
        ba_rows.append(dict(species=sp, n_pairs=len(br_list), bias_clr=bias, loa_lower=loa_l, loa_upper=loa_u, mean_abs_diff=madiff))
    write_tsv(OUT_SPECIES_BA, ba_rows, ["species","n_pairs","bias_clr","loa_lower","loa_upper","mean_abs_diff"])
    return OUT_SAMPLE, OUT_SPECIES_PREV, OUT_SPECIES_BA

if __name__ == '__main__':
    compute()
