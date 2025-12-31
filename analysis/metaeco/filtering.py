#!/usr/bin/env python3
# Portfolio-safe copy: identifiers generalized; runnable with included io_utils scaffold.

from __future__ import annotations
import math, json
from typing import List, Dict
from .io_utils import read_joined_table, write_tsv, float_or_nan, group_by, read_tsv

OUT_FLAGS = "results/comparisons/species_filter_flags.tsv"
OUT_CONC_FILT = "results/comparisons/species_concordance_filtered.tsv"
OUT_PREVBIN = "results/comparisons/species_concordance_prevbin_summary.tsv"
PARAMS_JSON = "results/comparisons/filtering_parameters_python.json"

PARAMS = dict(min_prevalence=2, min_max_fraction=1e-5, exclude_pattern_lower=["phytoplasma","phyllody","stunting"], drop_leading_apostrophe=True)

def compute():
    data = read_joined_table()
    by_species = group_by(data, 'species')
    flag_rows: List[Dict[str,object]] = []
    keep_species = []
    for sp, recs in by_species.items():
        prev_br = 0; prev_mp = 0; max_any = 0.0
        for r in recs:
            bf = float_or_nan(r.get('bracken_fraction'))
            ra = float_or_nan(r.get('rel_abundance'))
            if not math.isnan(bf) and bf>0: prev_br += 1; max_any = max(max_any, bf)
            if not math.isnan(ra) and ra>0: prev_mp += 1; max_any = max(max_any, ra)
        name_lower = sp.lower()
        flag_low_prev = max(prev_br, prev_mp) < PARAMS['min_prevalence']
        flag_low_abund = max_any < PARAMS['min_max_fraction']
        flag_pattern = any(pat in name_lower for pat in PARAMS['exclude_pattern_lower'])
        flag_leading_apostrophe = PARAMS['drop_leading_apostrophe'] and sp.startswith("'")
        exclude = flag_low_prev or flag_low_abund or flag_pattern or flag_leading_apostrophe
        if not exclude:
            keep_species.append(sp)
        flag_rows.append(dict(species=sp, prev_br=prev_br, prev_mp=prev_mp, max_any=max_any,
                              flag_low_prev=flag_low_prev, flag_low_abund=flag_low_abund,
                              flag_pattern=flag_pattern, flag_leading_apostrophe=flag_leading_apostrophe,
                              exclude=exclude))
    write_tsv(OUT_FLAGS, flag_rows, ["species","prev_br","prev_mp","max_any","flag_low_prev","flag_low_abund","flag_pattern","flag_leading_apostrophe","exclude"])
    # species-level concordance filtered
    by_sample = group_by(data, 'sampleID')
    from .io_utils import spearman
    spec_rows = []
    samples = sorted(by_sample.keys())
    sid_index = {sid:i for i,sid in enumerate(samples)}
    for sp in keep_species:
        recs = by_species[sp]
        br_vec = [0.0]*len(samples)
        mp_vec = [0.0]*len(samples)
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
        spec_rows.append(dict(species=sp, prev_br=prev_br, prev_mp=prev_mp, prev_any=prev_any, spearman_across_samples=rho))
    write_tsv(OUT_CONC_FILT, spec_rows, ["species","prev_br","prev_mp","prev_any","spearman_across_samples"])
    # prevalence bin summary
    bins = [(2,5),(6,10),(11,20),(21,30),(31,10**9)]
    prevbin_rows = []
    for lo,hi in bins:
        label = f"{lo}-{hi if hi<10**9 else '+'}" if hi<10**9 else f">={lo}"
        subset = [r for r in spec_rows if lo <= r['prev_any'] <= (hi if hi<10**9 else r['prev_any'])]
        if subset:
            vals = [r['spearman_across_samples'] for r in subset if not math.isnan(r['spearman_across_samples'])]
            if vals:
                import statistics
                prevbin_rows.append(dict(prev_bin=label, n_species=len(subset), median_rho=statistics.median(vals)))
    write_tsv(OUT_PREVBIN, prevbin_rows, ["prev_bin","n_species","median_rho"])
    # parameters
    with open(PARAMS_JSON,'w') as fh:
        import json; json.dump(PARAMS, fh, indent=2)
    return OUT_FLAGS, OUT_CONC_FILT, OUT_PREVBIN

if __name__ == '__main__':
    compute()
