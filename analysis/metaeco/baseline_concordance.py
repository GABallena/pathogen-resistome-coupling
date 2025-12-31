#!/usr/bin/env python3
# Portfolio-safe copy: identifiers generalized; runnable with included io_utils scaffold.

from __future__ import annotations
import math, os, json
from typing import Dict, List, Tuple
from .io_utils import read_joined_table, read_tsv, write_tsv, group_by, float_or_nan, clr_transform, geometric_mean, spearman

OUT_SAMPLE = "results/comparisons/sample_concordance_baseline.tsv"
OUT_SPECIES_PREV = "results/comparisons/species_prev_stratified_corr_baseline.tsv"
OUT_SPECIES_BA = "results/comparisons/species_bland_altman_clr_baseline.tsv"
OUT_QC = "results/comparisons/baseline_qc_summary.tsv"

CLR_VALUES = "results/comparisons/clr_values.tsv"

# Helper metrics

def lins_ccc(x: List[float], y: List[float]) -> float:
    import math
    if len(x) != len(y) or len(x) < 3:
        return math.nan
    xf = [v for v in x if not math.isnan(v)]
    yf = [v for v in y if not math.isnan(v)]
    if len(xf) < 3 or len(yf) < 3:
        return math.nan
    # align indices where both finite
    pair = [(a,b) for a,b in zip(x,y) if not(math.isnan(a) or math.isnan(b))]
    if len(pair) < 3: return math.nan
    xs = [a for a,_ in pair]; ys = [b for _,b in pair]
    mx = sum(xs)/len(xs); my = sum(ys)/len(ys)
    vx = sum((a-mx)**2 for a in xs)/(len(xs)-1)
    vy = sum((b-my)**2 for b in ys)/(len(ys)-1)
    sxy = sum((a-mx)*(b-my) for a,b in pair)/(len(xs)-1)
    denom = vx + vy + (mx - my)**2
    if denom == 0: return math.nan
    return (2 * sxy)/denom

def bland_altman(br: List[float], mp: List[float]) -> Tuple[float,float,float,float]:
    import math
    pair = [(a,b) for a,b in zip(br, mp) if not(math.isnan(a) or math.isnan(b))]
    if len(pair) < 3:
        return (math.nan, math.nan, math.nan, math.nan)
    diffs = [a-b for a,b in pair]
    bias = sum(diffs)/len(diffs)
    var = sum((d-bias)**2 for d in diffs)/(len(diffs)-1) if len(diffs)>1 else math.nan
    sd = math.sqrt(var) if not math.isnan(var) else math.nan
    loa_lower = bias - 1.96*sd if sd and not math.isnan(sd) else math.nan
    loa_upper = bias + 1.96*sd if sd and not math.isnan(sd) else math.nan
    return (bias, loa_lower, loa_upper, sum(abs(d) for d in diffs)/len(diffs))


def load_clr_dict():
    if not os.path.exists(CLR_VALUES):
        raise SystemExit("CLR values file missing; run clr_values module first.")
    rows = read_tsv(CLR_VALUES)
    # index by (sampleID, species)
    idx = {}
    for r in rows:
        sid = r.get('sampleID'); sp = r.get('species')
        cr = float_or_nan(r.get('clr_rel'))
        cb = float_or_nan(r.get('clr_bracken'))
        idx[(sid, sp)] = (cr, cb)
    return idx


def compute():
    joined = read_joined_table()
    clr_idx = load_clr_dict()
    by_sample = group_by(joined, 'sampleID')

    sample_rows: List[Dict[str,object]] = []
    for sid, recs in by_sample.items():
        # fractions & CLR alignment across overlapping species
        frac_pairs = []
        clr_pairs = []
        for r in recs:
            bf = float_or_nan(r.get('bracken_fraction'))
            ra = float_or_nan(r.get('rel_abundance'))
            sp = r['species']
            cr, cb = clr_idx.get((sid, sp), (math.nan, math.nan))
            frac_pairs.append((bf, ra))
            clr_pairs.append((cb, cr))  # order: bracken, metaphlan
        bf_list = [a for a,_ in frac_pairs]
        ra_list = [b for _,b in frac_pairs]
        cb_list = [a for a,_ in clr_pairs]
        cr_list = [b for _,b in clr_pairs]
        ccc_frac = lins_ccc(bf_list, ra_list)
        ccc_clr = lins_ccc(cb_list, cr_list)
        bias, loa_l, loa_u, madiff = bland_altman(cb_list, cr_list)
        sample_rows.append(dict(sampleID=sid, lins_ccc_fraction=ccc_frac, lins_ccc_clr=ccc_clr,
                                bland_altman_bias_clr=bias, bland_altman_loa_lower=loa_l, bland_altman_loa_upper=loa_u,
                                mean_abs_diff_clr=madiff))
    write_tsv(OUT_SAMPLE, sample_rows, ["sampleID","lins_ccc_fraction","lins_ccc_clr","bland_altman_bias_clr","bland_altman_loa_lower","bland_altman_loa_upper","mean_abs_diff_clr"])

    # prevalence stratified species spearman across samples
    by_species = group_by(joined, 'species')
    sample_list = sorted(by_sample.keys())
    sid_index = {sid:i for i,sid in enumerate(sample_list)}
    species_rows = []
    for sp, recs in by_species.items():
        br_vec = [0.0]*len(sample_list)
        mp_vec = [0.0]*len(sample_list)
        for r in recs:
            idx = sid_index[r['sampleID']]
            bf = float_or_nan(r.get('bracken_fraction'))
            ra = float_or_nan(r.get('rel_abundance'))
            if not math.isnan(bf): br_vec[idx] += bf
            if not math.isnan(ra): mp_vec[idx] += ra
        nz = [i for i,(a,b) in enumerate(zip(br_vec, mp_vec)) if a>0 or b>0]
        if len(nz) >= 3:
            rho = spearman([br_vec[i] for i in nz],[mp_vec[i] for i in nz])
        else:
            rho = math.nan
        prev_br = sum(1 for a in br_vec if a>0)
        prev_mp = sum(1 for a in mp_vec if a>0)
        prev_any = max(prev_br, prev_mp)
        # binning
        bins = [0,1,3,5,10,20,10**9]
        labels = ["0","1","2-3","4-5","6-10","11-20","21+"]
        for b1,b2,label in zip(bins[:-1], bins[1:], labels):
            if b1 <= prev_any < b2:
                prev_bin = label; break
        species_rows.append(dict(species=sp, spearman=rho, prev_br=prev_br, prev_mp=prev_mp, prev_bin=prev_bin))
    write_tsv(OUT_SPECIES_PREV, species_rows, ["species","spearman","prev_br","prev_mp","prev_bin"])

    # Species Bland-Altman CLR (intersection only per species across samples)
    clr_rows_input = read_tsv(CLR_VALUES)
    by_species_clr = group_by(clr_rows_input, 'species')
    ba_rows = []
    for sp, recs in by_species_clr.items():
        br_list=[]; mp_list=[]
        for r in recs:
            cb = float_or_nan(r.get('clr_bracken'))
            cr = float_or_nan(r.get('clr_rel'))
            if not math.isnan(cb) and not math.isnan(cr):
                br_list.append(cb); mp_list.append(cr)
        if len(br_list) < 1:
            continue
        bias, loa_l, loa_u, madiff = bland_altman(br_list, mp_list)
        ba_rows.append(dict(species=sp, n_pairs=len(br_list), bias_clr=bias, loa_lower=loa_l, loa_upper=loa_u, mean_abs_diff=madiff))
    write_tsv(OUT_SPECIES_BA, ba_rows, ["species","n_pairs","bias_clr","loa_lower","loa_upper","mean_abs_diff"])

    # QC summary
    qc = [dict(n_samples=len(by_sample), distinct_species_total=len(by_species), species_overlap_nonzero=sum(1 for s in by_species if s),
               baseline_species_prev_strata=len(species_rows), baseline_species_bland_altman=len(ba_rows))]
    write_tsv(OUT_QC, qc, ["n_samples","distinct_species_total","species_overlap_nonzero","baseline_species_prev_strata","baseline_species_bland_altman"])
    return OUT_SAMPLE, OUT_SPECIES_PREV, OUT_SPECIES_BA, OUT_QC

if __name__ == '__main__':
    compute()
