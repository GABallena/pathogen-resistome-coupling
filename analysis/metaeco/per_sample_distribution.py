#!/usr/bin/env python3
# Portfolio-safe copy: identifiers generalized; runnable with included io_utils scaffold.

"""
Per-sample distribution diagnostics for MetaPhlAn/Bracken outputs.

For each sample × rank × tool × transform, compute:
- n_taxa, n_nonzero, zero_fraction
- skewness (Fisher), excess kurtosis
- Jarque–Bera statistic and p-value (exact for chi-square df=2: p = exp(-JB/2))
- classification: parametric vs nonparametric (p>=alpha and zero_fraction<thresh_zero => parametric)
- suggested_correlation: Pearson (parametric) vs Spearman (nonparametric)

Transforms evaluated (auto-detected if columns exist):
- raw: abundance percent (tool column: metaphlan, bracken, bracken_taxonomic)
- clr: *_clr columns
- rclr: *_rclr columns
- logp: log10 of non-zero proportions (abundance/100)

Outputs a TSV summary.
"""
from __future__ import annotations

import argparse
import math
import os
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


def skew_kurtosis(x: np.ndarray) -> Tuple[float, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n < 3:
        return float("nan"), float("nan")
    mu = float(np.mean(x))
    v = x - mu
    m2 = float(np.mean(v ** 2))
    if m2 <= 0:
        return float("nan"), float("nan")
    m3 = float(np.mean(v ** 3))
    m4 = float(np.mean(v ** 4))
    g1 = m3 / (m2 ** 1.5)
    g2 = m4 / (m2 ** 2) - 3.0  # excess kurtosis
    return g1, g2


def jarque_bera(x: np.ndarray) -> Tuple[float, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n < 8:  # JB not meaningful on tiny n
        return float("nan"), float("nan")
    g1, g2 = skew_kurtosis(x)
    if not np.isfinite(g1) or not np.isfinite(g2):
        return float("nan"), float("nan")
    jb = n / 6.0 * (g1 ** 2 + (g2 ** 2) / 4.0)
    # For Chi-square with 2 DOF, survival function is exp(-x/2)
    p = math.exp(-jb / 2.0)
    return float(jb), float(p)


def build_summary(df: pd.DataFrame, rank: str, tool: str, col: str, transform: str, alpha: float, zero_thresh: float) -> pd.DataFrame:
    sub = df[df["rank"] == rank]
    if sub.empty or col not in sub.columns:
        return pd.DataFrame(columns=["sample","rank","tool","transform","n_taxa","n_nonzero","zero_fraction","skewness","excess_kurtosis","jb","pvalue","classification","suggested_correlation"])

    rows = []
    for sample, g in sub.groupby("sample"):
        vals = pd.to_numeric(g[col], errors="coerce")
        vals = vals[np.isfinite(vals)]
        n_taxa = int(vals.size)
        if n_taxa == 0:
            continue
        # Prepare series per transform
        if transform == "raw":
            x = vals.values
            zeros = np.sum(x <= 0)
        elif transform == "logp":
            # log10 of non-zero proportions
            p = (vals / 100.0).values
            mask = p > 0
            x = np.log10(p[mask])
            zeros = np.sum(~mask)
        elif transform in ("clr","rclr"):
            x = vals.values
            # Treat NaN as missing; zeros concept doesn't apply; approximate zeros as NaN
            zeros = int(np.sum(~np.isfinite(x)))
            x = x[np.isfinite(x)]
        else:
            continue
        n_nonzero = int(n_taxa - zeros)
        zero_fraction = float(zeros / n_taxa)
        g1, g2 = skew_kurtosis(x)
        jb, p = jarque_bera(x)
        is_param = (p >= alpha) and (zero_fraction < zero_thresh)
        classification = "parametric" if is_param else "nonparametric"
        suggested = "pearson" if is_param else "spearman"
        rows.append({
            "sample": sample,
            "rank": rank,
            "tool": tool,
            "transform": transform,
            "n_taxa": n_taxa,
            "n_nonzero": n_nonzero,
            "zero_fraction": zero_fraction,
            "skewness": g1,
            "excess_kurtosis": g2,
            "jb": jb,
            "pvalue": p,
            "classification": classification,
            "suggested_correlation": suggested,
        })
    return pd.DataFrame(rows)


def detect_tool_columns(df: pd.DataFrame) -> Dict[str, Dict[str, str]]:
    """Detect available columns per tool and transform.
    Returns dict tool -> transform -> column name.
    """
    mapping: Dict[str, Dict[str, str]] = {}
    # Base abundances
    tools = []
    if "metaphlan" in df.columns:
        tools.append("metaphlan")
    if "bracken" in df.columns:
        tools.append("bracken")
    if "bracken_taxonomic" in df.columns:
        tools.append("bracken_taxonomic")
    for t in tools:
        mapping[t] = {}
        mapping[t]["raw"] = t
        clr_col = f"{t}_clr" if t != "bracken_taxonomic" else None  # only have *_clr for original tools
        rclr_col = f"{t}_rclr" if t != "bracken_taxonomic" else None
        if clr_col and clr_col in df.columns:
            mapping[t]["clr"] = clr_col
        if rclr_col and rclr_col in df.columns:
            mapping[t]["rclr"] = rclr_col
        # logp will use the raw column name
        mapping[t]["logp"] = t
    return mapping


def main():
    ap = argparse.ArgumentParser(description="Per-sample distribution diagnostics (normality, skewness, kurtosis)")
    ap.add_argument("--input", default="diversity_results/metaphlan_bracken_merged_wide.tsv", help="Wide TSV with sample,taxon,rank, tool columns")
    ap.add_argument("--out", default="diversity_results/diagnostics/sample_distribution.tsv", help="Output TSV path")
    ap.add_argument("--ranks", default="all", help="Comma-separated ranks or 'all'")
    ap.add_argument("--transforms", default="raw,clr,rclr,logp", help="Comma-separated transforms to compute")
    ap.add_argument("--alpha", type=float, default=0.05, help="Significance threshold for normality (JB)")
    ap.add_argument("--zero-thresh", type=float, default=0.2, dest="zero_thresh", help="Zero fraction threshold; above this, classify nonparametric")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    df = pd.read_csv(args.input, sep="\t")
    if df.empty:
        print(f"ERROR: input is empty: {args.input}")
        return
    tools = detect_tool_columns(df)
    if not tools:
        print("ERROR: no tool columns detected in input.")
        return
    ranks = [r.strip() for r in (df["rank"].unique().tolist() if args.ranks == "all" else args.ranks.split(","))]
    transforms = [t.strip() for t in args.transforms.split(",") if t.strip()]

    out_frames: List[pd.DataFrame] = []
    for rk in ranks:
        for tool, txmap in tools.items():
            for tr in transforms:
                col = txmap.get(tr)
                if not col:
                    continue
                try:
                    out_frames.append(build_summary(df, rk, tool, col, tr, args.alpha, args.zero_thresh))
                except Exception as e:
                    print(f"WARN: failed summary for rank={rk}, tool={tool}, transform={tr}: {e}")
    out = pd.concat(out_frames, ignore_index=True) if out_frames else pd.DataFrame()
    if not out.empty:
        # Order columns
        cols = [
            "sample","rank","tool","transform","n_taxa","n_nonzero","zero_fraction",
            "skewness","excess_kurtosis","jb","pvalue","classification","suggested_correlation"
        ]
        out = out[cols]
        out.to_csv(args.out, sep="\t", index=False)
        print(f"Wrote: {args.out} (rows={len(out)})")
    else:
        print("WARN: nothing to write.")


if __name__ == "__main__":
    main()
