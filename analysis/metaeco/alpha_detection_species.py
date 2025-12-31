#!/usr/bin/env python3
# Portfolio-safe copy: identifiers generalized; runnable with included io_utils scaffold.

from __future__ import annotations
import math
from typing import Dict, List
from .io_utils import read_joined_table, group_by, float_or_nan, write_tsv

OUT_ALPHA = "results/comparisons/alpha_diversity.tsv"
OUT_DET = "results/comparisons/detection_overlap.tsv"
OUT_SPECIES_CONC = "results/comparisons/species_concordance.tsv"


def compute():
    rows = read_joined_table()
    by_sample = group_by(rows, 'sampleID')
    # Alpha diversity per tool
    alpha_rows: List[Dict[str,object]] = []
    # build per sample tool lists
    for sid, recs in by_sample.items():
        vals_br = [float_or_nan(r.get('bracken_fraction')) for r in recs]
        vals_ra = [float_or_nan(r.get('rel_abundance')) for r in recs]
        # species lists
        spp_br = [r['species'] for r,v in zip(recs, vals_br) if v and not math.isnan(v) and v>0]
        spp_mp = [r['species'] for r,v in zip(recs, vals_ra) if v and not math.isnan(v) and v>0]
        def alpha(vals, tool):
            vals = [v for v in vals if v and not math.isnan(v) and v>0]
            richness = len(vals)
            if richness>0:
                s = sum(vals)
                p = [v/s for v in vals]
                shannon = -sum(pi*math.log(pi) for pi in p if pi>0)
                simpson_d = sum(pi*pi for pi in p)
            else:
                shannon = math.nan; simpson_d = math.nan
            inv_simpson = (1/simpson_d) if simpson_d and not math.isnan(simpson_d) else math.nan
            pielou = (shannon/math.log(richness)) if richness>1 and not math.isnan(shannon) else math.nan
            def topN(th):
                if richness==0: return 0
                ordered = sorted(vals, reverse=True)
                cs = 0.0; total = sum(ordered)
                for i,v in enumerate(ordered):
                    cs += v
                    if cs/total >= th: return i+1
                return richness
            return dict(sampleID=sid, tool=tool, richness=richness, shannon=shannon, simpson_d=simpson_d,
                        inverse_simpson=inv_simpson, pielou=pielou, topN50=topN(0.5), topN80=topN(0.8), topN95=topN(0.95), total_abundance=sum(vals))
        alpha_rows.append(alpha(vals_br, 'Bracken'))
        alpha_rows.append(alpha(vals_ra, 'MetaPhlAn'))
    write_tsv(OUT_ALPHA, alpha_rows, ["sampleID","tool","richness","shannon","simpson_d","inverse_simpson","pielou","topN50","topN80","topN95","total_abundance"])

    # Detection overlap
    det_rows = []
    for sid, recs in by_sample.items():
        br = {r['species'] for r in recs if (v:=float_or_nan(r.get('bracken_fraction'))) and not math.isnan(v) and v>0}
        mp = {r['species'] for r in recs if (v:=float_or_nan(r.get('rel_abundance'))) and not math.isnan(v) and v>0}
        both = br & mp
        only_br = br - mp
        only_mp = mp - br
        union = len(both)+len(only_br)+len(only_mp)
        det_rows.append(dict(sampleID=sid, n_both=len(both), n_bracken_only=len(only_br), n_metaphlan_only=len(only_mp), union_species=union,
                              prop_both=(len(both)/union if union else math.nan), prop_bracken_only=(len(only_br)/union if union else math.nan), prop_metaphlan_only=(len(only_mp)/union if union else math.nan)))
    write_tsv(OUT_DET, det_rows, ["sampleID","n_both","n_bracken_only","n_metaphlan_only","union_species","prop_both","prop_bracken_only","prop_metaphlan_only"])

    # Species concordance (spearman across samples)
    # Build per species sample vectors
    by_species = group_by(rows, 'species')
    spec_rows = []
    from .io_utils import spearman
    # gather list of all sampleIDs for alignment
    samples = sorted(by_sample.keys())
    for sp, recs in by_species.items():
        v_br = [0.0]*len(samples)
        v_mp = [0.0]*len(samples)
        sid_index = {sid:i for i,sid in enumerate(samples)}
        for r in recs:
            idx = sid_index.get(r['sampleID'])
            if idx is None: continue
            bf = float_or_nan(r.get('bracken_fraction'))
            ra = float_or_nan(r.get('rel_abundance'))
            if not math.isnan(bf): v_br[idx] += bf
            if not math.isnan(ra): v_mp[idx] += ra
        nz = [i for i,(a,b) in enumerate(zip(v_br,v_mp)) if a>0 or b>0]
        if len(nz) >= 3:
            vx = [v_br[i] for i in nz]
            vy = [v_mp[i] for i in nz]
            rho = spearman(vx, vy)
        else:
            rho = math.nan
        spec_rows.append(dict(species=sp,
                              prevalence_bracken=sum(1 for a in v_br if a>0),
                              prevalence_metaphlan=sum(1 for a in v_mp if a>0),
                              mean_bracken=(sum(v_br)/len(v_br) if v_br else math.nan),
                              mean_metaphlan=(sum(v_mp)/len(v_mp) if v_mp else math.nan),
                              spearman_across_samples=rho))
    write_tsv(OUT_SPECIES_CONC, spec_rows, ["species","prevalence_bracken","prevalence_metaphlan","mean_bracken","mean_metaphlan","spearman_across_samples"])

    return OUT_ALPHA, OUT_DET, OUT_SPECIES_CONC

if __name__ == '__main__':
    compute()
