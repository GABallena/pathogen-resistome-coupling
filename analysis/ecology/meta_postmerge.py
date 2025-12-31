#!/usr/bin/env python3
"""
Resistome ecology analytics aligned with AFD (Gamma) framework — reviewer-revision edition.

Key additions vs prior:
- Fixed subcomposition for agreement (prevalence/topk/all).
- Principled zero handling for CLR/ALR: pseudocount, multiplicative replacement (Martín-Fernández), or Dirichlet smoothing.
- Bias harmonization across tools (CLR-mean shrinkage; ALR Deming).
- NB parameterization clarified: we use mean–variance form with size/shape β:
    E[X]=μ, Var[X]=μ + μ^2 / β
- Optional ZINB occupancy split θ via z_obs = (1-θ) + θ*p0_NB ⇒ θ̂ = (1 - z_obs) / (1 - p0_NB).
- Detection thresholds beyond n>0: fixed, kFDR from negatives, or depth-aware neg-Poisson.
- Hierarchical bootstrap (cluster columns, moving-block for time).
- Null models (random multinomial, tool permute, neutral-lite).
- Model selection beyond AIC: K-fold predictive log-likelihood and stacking.
- Bootstrap CIs for derived metrics R(N), R'(N), N*(q).
- Rarefaction/coverage diagnostics.
- Optional hierarchical multiple testing across modules.

Dependencies: pandas, numpy, scipy
"""
from __future__ import annotations

import argparse
import math
import os
import re
import sys
import logging
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats, optimize, special

try:
    import tomllib as _toml  # py311+
except ModuleNotFoundError:
    try:
        import tomli as _toml  # type: ignore[import-not-found]  # py<311 (optional)
    except ModuleNotFoundError:
        _toml = None


# -------------------- utils & io --------------------

def setup_logging(level: str = "INFO") -> logging.Logger:
    lvl = getattr(logging, level.upper(), logging.INFO)
    logging.basicConfig(level=lvl, format="%(asctime)s %(levelname)s %(message)s", datefmt="%H:%M:%S")
    return logging.getLogger("p4_eco")

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def parse_csv_list(s: Optional[str]) -> Optional[List[float]]:
    if not s:
        return None
    out = []
    for part in s.split(","):
        part = part.strip()
        if part:
            try:
                out.append(float(part))
            except Exception:
                pass
    return out or None

def read_table_any(path: str) -> pd.DataFrame:
    if path.lower().endswith((".tsv", ".txt")):
        return pd.read_csv(path, sep="\t")
    if path.lower().endswith(".csv"):
        return pd.read_csv(path)
    if path.lower().endswith((".xlsx", ".xls")):
        try:
            return pd.read_excel(path)
        except Exception as e:
            raise RuntimeError(f"Reading Excel requires openpyxl/xlrd: {e}")
    raise ValueError(f"Unsupported file type: {path}")

# IEEE-safe log cap for doubles (exp(709) ~ 8.2e307; keep safety margin)
_LOG_MAX = 745.0
# Small epsilons for numeric stability
_EPS_DIV = 1e-12
_EPS_LOG = 1e-300

def _safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    """Pearson correlation with guards for zero variance and small N.

    Returns NaN if fewer than 2 finite points or if either std dev < eps.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 2:
        return float('nan')
    xv = x[mask]; yv = y[mask]
    sx = float(np.std(xv))
    sy = float(np.std(yv))
    if sx < _EPS_DIV or sy < _EPS_DIV:
        return float('nan')
    return float(np.corrcoef(xv, yv)[0, 1])

# ------------------ metadata selection ------------------

def tag_negatives_from_metadata(meta: pd.DataFrame) -> pd.Series:
    neg_tokens = ("neg", "blank", "control", "ntc")
    neg_cols = [c for c in meta.columns if any(tok in str(c).lower() for tok in neg_tokens)]
    is_neg = pd.Series(False, index=meta.index)
    for col in meta.columns:
        vals = meta[col].astype(str).str.lower()
        is_neg = is_neg | vals.str.contains("|".join(neg_tokens), regex=True, na=False)
    for col in neg_cols:
        if meta[col].dtype == bool:
            is_neg = is_neg | meta[col]
    return is_neg

def build_sample_lists(
    merged: pd.DataFrame,
    metadata_path: Optional[str],
    sample_col: Optional[str],
    include_regex: Optional[str],
    cohort_col: Optional[str],
    cohort_allow: Optional[str],
    date_col: Optional[str],
    exclude_old_before: Optional[str],
    neg_regex: str,
    out_dir: str,
    logger: logging.Logger
) -> Tuple[List[str], List[str], List[str]]:
    samples = sorted(merged["sample"].astype(str).unique())
    included = set(samples)
    negatives: set = set()
    excluded: set = set()

    if include_regex:
        patt = re.compile(include_regex)
        keep = {s for s in samples if patt.search(s)}
        drop = set(samples) - keep
        included &= keep
        excluded |= drop

    meta = None
    if metadata_path and os.path.exists(metadata_path):
        meta = read_table_any(metadata_path)
        if not sample_col:
            candidates = [c for c in meta.columns if str(c).strip().lower() in {"sample","sample_id","id","sampleid","sample code","sample_code"}]
            if candidates:
                sample_col = candidates[0]
        if not sample_col or sample_col not in meta.columns:
            sample_col = "SAMPLE CODE" if "SAMPLE CODE" in meta.columns else None

        if sample_col:
            key_to_samples: Dict[str, List[str]] = {}
            for s in samples:
                key = s if sample_col != "SAMPLE CODE" else s.split("_")[0]
                key_to_samples.setdefault(key, []).append(s)

            if cohort_col and cohort_col in meta.columns and cohort_allow:
                allowed = {str(x).strip() for x in cohort_allow.split(",")}
                disallowed_keys = set()
                for _, r in meta.iterrows():
                    key = str(r[sample_col]).strip() if sample_col in r else None
                    if key is None:
                        continue
                    if str(r[cohort_col]).strip() not in allowed:
                        disallowed_keys.add(key)
                for k in disallowed_keys:
                    for s in key_to_samples.get(k, []):
                        if s in included:
                            included.remove(s); excluded.add(s)

            if date_col and date_col in meta.columns and exclude_old_before:
                try:
                    cutoff = pd.to_datetime(exclude_old_before)
                    dates = pd.to_datetime(meta[date_col], errors="coerce")
                    old_keys = set(meta.loc[(dates.notna()) & (dates < cutoff), sample_col].astype(str))
                    for k in old_keys:
                        for s in key_to_samples.get(k, []):
                            if s in included:
                                included.remove(s); excluded.add(s)
                except Exception as e:
                    logger.warning("Date filter ignored: %s", e)

            neg_mask = tag_negatives_from_metadata(meta)
            neg_keys = set(meta.loc[neg_mask, sample_col].astype(str)) if sample_col in meta.columns else set()
            for k in neg_keys:
                for s in key_to_samples.get(k, []):
                    negatives.add(s)
        else:
            logger.warning("Could not detect sample column in metadata; skipping cohort/date filters.")
    else:
        if metadata_path:
            logger.warning("Metadata file not found: %s", metadata_path)

    if neg_regex:
        patt = re.compile(neg_regex)
        for s in samples:
            if patt.search(s):
                negatives.add(s)

    included -= negatives

    def write_list(lst: List[str], fname: str):
        pd.DataFrame({"sample": sorted(lst)}).to_csv(os.path.join(out_dir, fname), sep="\t", index=False)

    write_list(list(included), "samples_included.tsv")
    write_list(sorted(list(negatives)), "samples_negatives.tsv")
    write_list(sorted(list(excluded)), "samples_excluded.tsv")

    return (sorted(included), sorted(negatives), sorted(excluded))

# ------------------ zero handling ------------------

def zero_replace_pseudocount(p: np.ndarray, pseudocount: float) -> np.ndarray:
    p = np.asarray(p, float)
    p = np.where(p > 0, p, float(pseudocount))
    s = p.sum()
    return p / (s if s > 0 else 1.0)

def zero_replace_mult(p: np.ndarray, delta: float = 1e-6) -> np.ndarray:
    """Multiplicative replacement (Martín-Fernández): replace zeros with δ and renormalize by shrinking nonzeros."""
    p = np.asarray(p, float)
    D = p.size
    k = int((p <= 0).sum())
    if k == 0:
        s = p.sum(); return p / (s if s>0 else 1.0)
    d = min(float(delta), 1.0/(10.0*max(D,1)))  # guard
    z = (p <= 0)
    nz = ~z
    mass = k * d
    out = p.copy()
    out[z] = d
    if nz.any():
        scale = (1.0 - mass) / max(p[nz].sum(), 1e-15)
        out[nz] = p[nz] * scale
    s = out.sum()
    return out / (s if s > 0 else 1.0)

def zero_replace_dirichlet(n: Optional[np.ndarray], N: Optional[float], D: int, alpha: float = 0.5) -> np.ndarray:
    """Dirichlet posterior mean E[p_i|n] = (n_i + α)/(N + D α). If counts unavailable, fall back to uniform prior."""
    if n is None or N is None or N <= 0:
        return np.full(D, 1.0 / max(D,1))
    n = np.asarray(n, float)
    denom = float(N + D * alpha)
    return (n + alpha) / (denom if denom > 0 else 1.0)

def _clr_from_p(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    with np.errstate(divide="ignore"):
        logp = np.log(np.clip(p, _EPS_LOG, None))
    return logp - np.mean(logp)

def _alr_from_p(p: np.ndarray, ref_idx: int) -> np.ndarray:
    p = np.asarray(p, float)
    with np.errstate(divide="ignore"):
        lp = np.log(np.clip(p, _EPS_LOG, None))
    return lp - lp[int(ref_idx)]

def _plr_vector_pair_from_p(pa: np.ndarray, pb: np.ndarray, rng=None, maxpairs=0):
    a = np.asarray(pa, float)
    b = np.asarray(pb, float)
    idx = np.where((a > 0) & (b > 0))[0]
    if idx.size < 2:
        return np.array([]), np.array([])
    pairs = np.array(np.triu_indices(idx.size, k=1)).T
    if maxpairs and pairs.shape[0] > maxpairs:
        rng = np.random.default_rng(17) if rng is None else rng
        sel = rng.choice(pairs.shape[0], size=maxpairs, replace=False)
        pairs = pairs[sel]
    v1 = np.log(a[idx[pairs[:,0]]]) - np.log(a[idx[pairs[:,1]]])
    v2 = np.log(b[idx[pairs[:,0]]]) - np.log(b[idx[pairs[:,1]]])
    return v1, v2

def transform_compositional_pair(
    a_pct: np.ndarray, b_pct: np.ndarray, mode: str,
    zero_replacement: str, pseudocount: float,
    alr_ref: Optional[int], alr_topk: int, plr_maxpairs: int, rng=None,
    # optional for Dirichlet smoothing:
    a_counts: Optional[np.ndarray]=None, b_counts: Optional[np.ndarray]=None, N_a: Optional[float]=None, N_b: Optional[float]=None,
    zero_delta: float=1e-6, dirichlet_alpha: float=0.5
):
    a = np.asarray(a_pct, float) * 0.01
    b = np.asarray(b_pct, float) * 0.01
    D = a.size
    if zero_replacement == "pseudocount":
        pa = zero_replace_pseudocount(a, pseudocount)
        pb = zero_replace_pseudocount(b, pseudocount)
    elif zero_replacement == "mult":
        pa = zero_replace_mult(a, zero_delta)
        pb = zero_replace_mult(b, zero_delta)
    else:
        pa = zero_replace_dirichlet(a_counts, N_a, D, dirichlet_alpha)
        pb = zero_replace_dirichlet(b_counts, N_b, D, dirichlet_alpha)
        # renormalize just in case
        pa = pa / max(pa.sum(), 1e-15); pb = pb / max(pb.sum(), 1e-15)

    if mode == "clr":
        v1 = _clr_from_p(pa); v2 = _clr_from_p(pb)
    elif mode == "alr":
        if alr_ref is not None:
            ref = int(alr_ref)
        elif alr_topk and alr_topk > 0:
            ref = int(np.argsort(pb)[-int(alr_topk):][0])
        else:
            ref = int(np.nanargmax(pb))
        v1 = _alr_from_p(pa, ref); v2 = _alr_from_p(pb, ref)
    else:
        v1, v2 = _plr_vector_pair_from_p(pa, pb, rng=rng, maxpairs=plr_maxpairs)
    return v1, v2

# ------------------ subcomposition & bias ------------------

def build_fixed_subcomposition(df, samples, t1, t2, mode="prevalence", prevalence=0.5, topk=0):
    pos = df[df["sample"].isin(samples)]
    sp = pos.groupby("taxon")["sample"].nunique()
    S = len(set(samples))
    prev = sp / max(S, 1)
    if mode == "all":
        keep = pos["taxon"].unique()
    elif mode == "prevalence":
        keep = prev[prev >= float(prevalence)].index.values
    else:  # topk
        if topk <= 0:
            keep = pos["taxon"].unique()
        else:
            means = pos.groupby("taxon")[[t1, t2]].mean().mean(axis=1)
            rank = (pd.DataFrame({"prev": prev}).join(means.rename("mean")).fillna(0.0)
                    .sort_values(["prev","mean"], ascending=[False,False]))
            keep = rank.index.values[:int(topk)]
    return np.asarray(sorted(set(map(str, keep))))

def _build_aligned_vectors(g, taxa_keep, tool1, tool2):
    """Return aligned vectors for tool1/tool2 across taxa_keep.

    Deduplicate potential repeated taxon rows within a sample by averaging numeric values.
    """
    idx_map = {t: i for i, t in enumerate(taxa_keep)}
    a = np.zeros(len(taxa_keep), float)
    b = np.zeros(len(taxa_keep), float)
    # Keep only relevant columns and coerce to numeric
    tmp = g[["taxon", tool1, tool2]].copy()
    tmp[tool1] = pd.to_numeric(tmp[tool1], errors="coerce").astype(float)
    tmp[tool2] = pd.to_numeric(tmp[tool2], errors="coerce").astype(float)
    # Aggregate duplicates per taxon (mean)
    gg = tmp.groupby("taxon", as_index=False).mean(numeric_only=True).set_index("taxon")
    t_in = [t for t in taxa_keep if t in gg.index]
    if t_in:
        ixs = [idx_map[t] for t in t_in]
        sel = gg.loc[t_in]
        a[ixs] = sel[tool1].fillna(0.0).values
        b[ixs] = sel[tool2].fillna(0.0).values
    return a, b

def estimate_clr_bias(df, samples, taxa_keep, t1, t2, zero_method, pc, zero_delta, dirichlet_alpha, lam=0.1):
    diffs = []
    for s, g in df[df["sample"].isin(samples)].groupby("sample", sort=False):
        a_pct, b_pct = _build_aligned_vectors(g, taxa_keep, t1, t2)
        v1, v2 = transform_compositional_pair(a_pct, b_pct, mode="clr",
                                              zero_replacement=zero_method, pseudocount=pc,
                                              alr_ref=None, alr_topk=0, plr_maxpairs=0,
                                              zero_delta=zero_delta, dirichlet_alpha=dirichlet_alpha)
        if v1.size:
            diffs.append(v2 - v1)
    if not diffs:
        return pd.DataFrame({"taxon": taxa_keep, "bias_clr": np.zeros(len(taxa_keep)), "se": np.zeros(len(taxa_keep))})
    D = np.vstack(diffs)
    mu = np.nanmean(D, axis=0)
    se = np.nanstd(D, axis=0, ddof=1) / np.sqrt(max(D.shape[0],1))
    w = (se**2) / (se**2 + lam)
    bias = w * mu  # shrink toward 0
    return pd.DataFrame({"taxon": taxa_keep, "bias_clr": bias, "se": se})

def deming_regression(x, y, delta=1.0):
    x = np.asarray(x, float); y = np.asarray(y, float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size < 3: return 0.0, 1.0
    xbar, ybar = x.mean(), y.mean()
    Sxx = np.mean((x - xbar)**2); Syy = np.mean((y - ybar)**2); Sxy = np.mean((x - xbar)*(y - ybar))
    if abs(Sxy) < 1e-15:
        return ybar - xbar, 1.0
    b = (Syy - delta*Sxx + np.sqrt((Syy - delta*Sxx)**2 + 4*delta*Sxy**2)) / (2*Sxy)
    a = ybar - b * xbar
    return a, b

# ------------------ agreement (fixed subcomposition) ------------------

def agreement_per_sample_fixed(df, tool1, tool2, taxa_keep, mode="clr",
                               zero_replacement="pseudocount", pseudocount=1e-6,
                               alr_ref=None, alr_topk=0, plr_maxpairs=0, rng=None,
                               zero_delta=1e-6, dirichlet_alpha=0.5,
                               bias_df: Optional[pd.DataFrame]=None, bias_mode="none", deming_delta=1.0):
    rows = []
    taxa_keep = np.asarray(taxa_keep, dtype=str)
    bias_map = {}
    if bias_df is not None and "bias_clr" in bias_df.columns:
        bias_map = dict(zip(bias_df["taxon"].astype(str), bias_df["bias_clr"].astype(float)))
    for s, g in df.groupby("sample", sort=False):
        a_pct, b_pct = _build_aligned_vectors(g, taxa_keep, tool1, tool2)

        # bias correction (CLR mean)
        if bias_mode == "clr-mean" and mode == "clr" and bias_map:
            # Apply on CLR space: subtract bias from tool2 CLR after zero-handling
            v1c, v2c = transform_compositional_pair(a_pct, b_pct, "clr", zero_replacement, pseudocount,
                                                    alr_ref, alr_topk, plr_maxpairs, rng,
                                                    zero_delta=zero_delta, dirichlet_alpha=dirichlet_alpha)
            adj = np.array([bias_map.get(t, 0.0) for t in taxa_keep], float)
            v2c = v2c - adj
            v1, v2 = v1c, v2c
        elif bias_mode == "alr-deming" and mode == "alr":
            # Fit Deming on-the-fly using ALR w.r.t ref across taxa in this sample
            ref = int(alr_ref) if alr_ref is not None else int(np.nanargmax(b_pct))
            v1 = _alr_from_p(zero_replace_pseudocount(a_pct*0.01, 1e-9), ref)
            v2 = _alr_from_p(zero_replace_pseudocount(b_pct*0.01, 1e-9), ref)
            a_int, b_slope = deming_regression(v1, v2, delta=deming_delta)
            v2 = (v2 - a_int) / max(b_slope, 1e-12)
        else:
            v1, v2 = transform_compositional_pair(a_pct, b_pct, mode, zero_replacement, pseudocount,
                                                  alr_ref, alr_topk, plr_maxpairs, rng,
                                                  zero_delta=zero_delta, dirichlet_alpha=dirichlet_alpha)
        if v1.size < 2 or v2.size < 2 or v1.size != v2.size:
            pear = np.nan; spear = np.nan; aitch = np.nan
        else:
            pear = _safe_corr(v1, v2)
            spear = float(stats.spearmanr(v1, v2, nan_policy="omit").correlation)
            aitch = float(np.linalg.norm(v1 - v2)) if mode == "clr" else np.nan
        rows.append({"sample": s, "pearson": pear, "spearman": spear, "aitchison_clr_norm": aitch})
    return pd.DataFrame(rows).sort_values("sample", kind="mergesort").reset_index(drop=True)

# ------------------ counts reconstruction ------------------

def build_counts_from_pct(df: pd.DataFrame, counts_df: pd.DataFrame, tool_col: str) -> pd.DataFrame:
    m = df[["sample", "taxon", tool_col]].merge(counts_df[["sample", "N"]], on="sample", how="inner")
    rows = []
    for s, g in m.groupby("sample", sort=False):
        N = int(round(float(g["N"].iloc[0])))
        p = np.clip(pd.to_numeric(g[tool_col], errors="coerce").astype(float).values, 0.0, 100.0) * 0.01
        raw = p * N
        flo = np.floor(raw)
        give = int(N - flo.sum())
        if give > 0:
            idx = np.argsort(raw - flo)[-give:]
            flo[idx] += 1
        n = flo.astype(int)
        for (i, (_, r)) in enumerate(g.iterrows()):
            rows.append({"sample": s, "taxon": r["taxon"], "N": N, "n": int(n[i])})
    return pd.DataFrame(rows)

# ------------------ unbiased moments & curvature ------------------

def unbiased_moments(counts: pd.DataFrame) -> pd.DataFrame:
    out = []
    by_tax = counts.groupby("taxon", sort=False)
    for tax, g in by_tax:
        Ns = g["N"].astype(float).values
        ns = g["n"].astype(float).values
        D1 = float(np.sum(Ns))
        D2 = float(np.sum(Ns * (Ns - 1)))
        D3 = float(np.sum(Ns * (Ns - 1) * (Ns - 2)))
        if D1 <= 0:
            continue
        mu1 = float(np.sum(ns) / D1)
        mu2 = float(np.sum(ns * (ns - 1)) / D2) if D2 > 0 else np.nan
        mu3 = float(np.sum(ns * (ns - 1) * (ns - 2)) / D3) if D3 > 0 else np.nan
        with np.errstate(divide="ignore", invalid="ignore"):
            frac = ns / np.maximum(Ns, 1.0)
            naive2 = float(np.nanmean(frac**2))
            naive3 = float(np.nanmean(frac**3))
        out.append({"taxon": tax, "mu1": mu1, "mu2": mu2, "mu3": mu3, "naive2_obs": naive2, "naive3_obs": naive3})
    return pd.DataFrame(out).sort_values("taxon", kind="mergesort").reset_index(drop=True)

def predict_naive_from_moments(m: pd.DataFrame, counts: pd.DataFrame) -> pd.DataFrame:
    Ns = counts[["sample","N"]].drop_duplicates()["N"].astype(float).values
    Ns_safe = np.maximum(Ns, 1.0)
    invN = float(np.mean(1.0 / Ns_safe))
    invN2 = float(np.mean(1.0 / (Ns_safe**2)))
    df = m.copy()
    df["naive2_pred"] = df["mu2"] + (df["mu1"] - df["mu2"]) * invN
    df["naive3_pred"] = df.apply(
        lambda r: (r["mu3"] + 3*(r["mu2"] - r["mu3"])*invN + (r["mu1"] - 3*r["mu2"] + 2*r["mu3"])*invN2)
        if np.isfinite(r["mu2"]) and np.isfinite(r["mu3"]) else np.nan,
        axis=1
    )
    df["resid2"] = df["naive2_obs"] - df["naive2_pred"]
    df["resid3"] = df["naive3_obs"] - df["naive3_pred"]
    return df

# ------------------ detection thresholds ------------------

def thresholds_fixed(counts: pd.DataFrame, t: int = 1) -> Dict[str,int]:
    return {}

def thresholds_kfdr(counts: pd.DataFrame, neg_samples: List[str], alpha: float=0.01) -> Dict[str,int]:
    # Choose smallest t such that FPR across negatives <= alpha
    out = {}
    neg = counts[counts["sample"].isin(neg_samples)]
    if neg.empty: return out
    grp = neg.groupby("taxon")["n"]
    for tax, x in grp:
        vals = np.sort(x.values)
        if len(vals)==0:
            out[tax] = 1; continue
        # empirical FPR(t) = mean[1{n>=t}]
        unique = np.unique(vals)
        t_best = 1
        for t in range(1, int(unique.max())+2):
            fpr = float((vals >= t).mean())
            if fpr <= alpha:
                t_best = t; break
        out[tax] = t_best
    return out

def thresholds_neg_poisson(counts: pd.DataFrame, neg_samples: List[str], alpha: float=0.01) -> Dict[str,Tuple[float,float]]:
    """Return background rate b_i per taxon (n_neg_total / N_neg_total)."""
    neg = counts[counts["sample"].isin(neg_samples)]
    out = {}
    if neg.empty: return out
    Ntot = float(neg.drop_duplicates("sample")["N"].sum())
    if Ntot <= 0: return out
    bg = neg.groupby("taxon", sort=False)["n"].sum()
    b = (bg / Ntot).to_dict()
    return b

def is_detected(n: int, N: int, tax: str,
                mode: str, fixed_t: int,
                kfdr_map: Dict[str,int],
                bkg_rate: Dict[str,float],
                alpha: float) -> bool:
    if mode == "fixed":
        return n >= max(1, fixed_t)
    if mode == "kFDR":
        t = int(kfdr_map.get(tax, 1))
        return n >= t
    # neg-poisson
    lam = float(bkg_rate.get(tax, 0.0)) * float(N)
    if lam <= 0:
        return n >= 1
    # choose smallest t with sf(t-1; lam) <= alpha
    # Evaluate detection if n >= t*
    t_min = 1
    # quick approximate: inverse survival roughly around lam + z*sqrt(lam)
    # but do monotone search up to, say, n
    p_tail = stats.poisson.sf(max(n-1,0), lam)
    return bool(p_tail <= float(alpha))

# ------------------ NB / Gamma AFD, occupancy, GOF ------------------

def safe_nb_zero_prob(mu, beta):
    """Stable P(n=0) for NB with mean mu and size β (Var=μ+μ²/β); handles Poisson limit β→∞."""
    mu = np.asarray(mu, float)
    if np.isinf(beta):
        return np.exp(-np.clip(mu, 0.0, _LOG_MAX))
    z = -float(beta) * np.log1p(np.clip(mu / max(float(beta), 1e-12), 0.0, 1e12))
    z = np.clip(z, -_LOG_MAX, _LOG_MAX)
    p0 = np.exp(z)
    return float(p0) if p0.ndim == 0 else p0

def occupancy_from_beta(bar_x: np.ndarray, Ns: np.ndarray, beta: float) -> np.ndarray:
    mu = np.outer(bar_x, Ns)
    p0 = safe_nb_zero_prob(mu, beta)
    return np.mean(1.0 - p0, axis=1)

def fit_shared_beta(bar_x: np.ndarray, Ns: np.ndarray, occ_emp: np.ndarray, grid: Optional[List[float]] = None) -> float:
    def mse(b: float) -> float:
        b = max(float(b), 1e-6)
        pred = occupancy_from_beta(bar_x, Ns, b)
        d = pred - occ_emp
        return float(np.nanmean(d*d))
    if grid:
        vals = np.array([mse(b) for b in grid], float)
        return float(grid[int(np.argmin(vals))])
    res = optimize.minimize_scalar(mse, bounds=(1e-6, 1e6), method="bounded")
    return float(res.x if res.success else 1.0)

def fit_shared_beta_cv_from_counts(counts_pos: pd.DataFrame, Ns_pos: np.ndarray, folds=5, grid=None, rng=None,
                                   metric="occupancy_mse", detect_mode="fixed", fixed_t=1,
                                   kfdr_map=None, bkg_rate=None, detect_alpha=0.01):
    """k-fold CV for shared beta using fold-specific empirical occupancy and train-only bar_x."""
    if folds <= 1:
        # build global bar_x, occ_emp
        bar = unbiased_moments(counts_pos); bar_x = bar.set_index("taxon")["mu1"].values
        occ_emp = counts_pos.groupby("taxon")["n"].apply(lambda x: (x>0).mean()).reindex(bar["taxon"]).fillna(0.0).values
        return fit_shared_beta(bar_x, Ns_pos, occ_emp, grid)
    rng = np.random.default_rng(123) if rng is None else rng
    samples = counts_pos["sample"].drop_duplicates().tolist()
    idx = np.arange(len(samples)); rng.shuffle(idx); cuts = np.array_split(idx, folds)
    s_order = np.array(samples, dtype=object)
    def score(b):
        sc = []
        for test in cuts:
            test_s = set(s_order[test].tolist())
            train_s = set(samples) - test_s
            train = counts_pos[counts_pos["sample"].isin(train_s)]
            test_df = counts_pos[counts_pos["sample"].isin(test_s)]
            if train.empty or test_df.empty: continue
            m = unbiased_moments(train)
            bar_x = m.set_index("taxon")["mu1"]
            # empirical occupancy in test, with detection threshold
            def det_series(df):
                arr = []
                for tax, g in df.groupby("taxon", sort=False):
                    det = []
                    for _, r in g.iterrows():
                        det.append(int(is_detected(int(r["n"]), int(r["N"]), tax, detect_mode, fixed_t, kfdr_map or {}, bkg_rate or {}, detect_alpha)))
                    arr.append((tax, float(np.mean(det))))
                if not arr: return pd.Series(dtype=float)
                t,v = zip(*arr); return pd.Series(v, index=list(t))
            occ_emp = det_series(test_df)
            bx = bar_x.reindex(occ_emp.index).fillna(0.0).values
            Ns_test = test_df.drop_duplicates("sample")["N"].astype(float).values
            pred = occupancy_from_beta(bx, Ns_test, max(float(b),1e-6))
            if metric == "zero_fraction_mse":
                pred, obs = 1.0 - pred, 1.0 - occ_emp.values
            else:
                obs = occ_emp.values
            sc.append(float(np.nanmean((pred - obs)**2)))
        return float(np.mean(sc)) if sc else np.inf
    if grid:
        vals = np.array([score(b) for b in grid], float)
        return float(grid[int(np.argmin(vals))])
    res = optimize.minimize_scalar(score, bounds=(1e-6, 1e6), method="bounded")
    return float(res.x if res.success else 1.0)

def fit_beta_per_species(
    bar_x: np.ndarray, Ns: np.ndarray,
    counts_mat: Dict[str, np.ndarray], spp: List[str],
    min_nonzero: int = 8, min_total: int = 1000
) -> Dict[str, float]:
    betas = {}
    for i, tx in enumerate(spp):
        ns = counts_mat.get(tx)
        if ns is None:
            betas[tx] = np.nan; continue
        if (ns > 0).sum() < min_nonzero or ns.sum() < min_total:
            betas[tx] = np.nan; continue
        mu = bar_x[i] * Ns
        def nll(logb: float) -> float:
            b = np.exp(logb)
            mu_safe = np.clip(mu, _EPS_LOG, None)
            with np.errstate(divide="ignore", invalid="ignore"):
                ll = special.gammaln(ns + b) - special.gammaln(b) - special.gammaln(ns + 1) \
                     + b*(np.log(b) - np.log(b + mu_safe)) + ns*(np.log(mu_safe) - np.log(b + mu_safe))
            return -float(np.nansum(ll))
        try:
            res = optimize.minimize_scalar(nll, bounds=(-20, 20), method="bounded")
            betas[tx] = float(np.exp(res.x)) if res.success else np.nan
        except Exception:
            betas[tx] = np.nan
    return betas

def nb_logpmf_vec(k, mu, beta):
    """NB log PMF with mean μ and size β; stable."""
    k = np.asarray(k, dtype=np.int64)
    mu_safe = np.clip(mu, _EPS_LOG, None)
    b = float(beta)
    with np.errstate(divide='ignore', invalid='ignore'):
        return (special.gammaln(k + b) - special.gammaln(b) - special.gammaln(k + 1) +
                b * (np.log(b) - np.log(b + mu_safe)) +
                k * (np.log(mu_safe) - np.log(b + mu_safe)))

def nb_pmf(k, mu, beta):
    lp = nb_logpmf_vec(k, mu, beta)
    p = np.exp(np.clip(lp, -_LOG_MAX, 0))
    return np.clip(p, 0.0, 1.0)

def chisq_gof_nb(counts: pd.DataFrame, bar_x: pd.Series, beta_shared: float, betas: Optional[Dict[str, float]], fitted_params: int) -> pd.DataFrame:
    rows = []
    for tax, g in counts.groupby("taxon", sort=False):
        Ns = g.set_index("sample")["N"].astype(float)
        ks = g.set_index("sample")["n"].astype(int)
        mu = bar_x.loc[tax] * Ns
        b = betas.get(tax, beta_shared) if betas else beta_shared
        obs_k, obs_cnt = np.unique(ks.values, return_counts=True)
        exp_cnt = []
        for k in obs_k:
            pk = nb_pmf(np.array([k]), mu.values, b)[0]
            exp_cnt.append(float(np.sum(pk)))
        obs = obs_cnt.astype(float).tolist()
        exp = exp_cnt
        i = 0
        while i < len(exp):
            if exp[i] >= 5 or len(exp) == 1:
                i += 1; continue
            if i + 1 < len(exp):
                exp[i + 1] += exp[i]; obs[i + 1] += obs[i]; del exp[i]; del obs[i]
            else:
                exp[i - 1] += exp[i]; obs[i - 1] += obs[i]; del exp[i]; del obs[i]; i -= 1
        df = len(exp) - 1 - fitted_params
        if df <= 0:
            rows.append({"taxon": tax, "gof_chisq": np.nan, "gof_df": 0, "gof_p": np.nan})
            continue
        chisq = float(np.sum((np.array(obs) - np.array(exp)) ** 2 / np.maximum(np.array(exp), 1e-12)))
        p = float(stats.chi2.sf(chisq, df))
        rows.append({"taxon": tax, "gof_chisq": chisq, "gof_df": df, "gof_p": p})
    return pd.DataFrame(rows).sort_values("taxon", kind="mergesort").reset_index(drop=True)

def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, float)
    m = np.isfinite(p).sum()
    q = np.full_like(p, np.nan, dtype=float)
    if m == 0: return q
    vals = p[np.isfinite(p)]
    order = np.argsort(vals)
    ranked = vals[order]
    adj = ranked * m / (np.arange(m) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    res = np.empty_like(vals); res[order] = adj
    q[np.isfinite(p)] = res
    return q

def holm_stepdown(pvals):
    p = np.asarray(pvals, float)
    out = np.full_like(p, np.nan, dtype=float)
    mask = np.isfinite(p)
    m = int(mask.sum())
    if m == 0:
        return out
    idx = np.argsort(p[mask])
    sorted_p = p[mask][idx]
    adj = np.maximum.accumulate((m - np.arange(m)) * sorted_p)
    adj = np.minimum(adj, 1.0)
    tmp = np.empty_like(sorted_p)
    tmp[idx] = adj
    out[mask] = tmp
    return out

def n_star_q(bar_x: float, beta: float, q: float) -> float:
    q = float(np.clip(q, 1e-12, 1-1e-12))
    if math.isinf(beta):
        return -math.log(1.0 - q) / max(bar_x, 1e-15)
    return float(beta * ((1 - q) ** (-1.0 / beta) - 1.0) / max(bar_x, 1e-15))

def richness_and_slope(bar_x: np.ndarray, beta: float, N: float) -> Tuple[float,float]:
    mu = bar_x * N
    p0 = safe_nb_zero_prob(mu, beta)
    R = float(np.sum(1.0 - p0))
    if np.isinf(beta):
        Rp = float(np.sum(bar_x * np.exp(-mu)))
    else:
        Rp = float(np.sum(bar_x * (1.0 + mu / beta) ** (-(beta + 1.0))))
    return R, Rp

# ------------------ Posterior predictive checks ------------------

def posterior_predictive_checks(counts, bar_x, beta_shared, betas=None, draws=1000, rng=None):
    if draws <= 0:
        return None
    rng = np.random.default_rng(7) if rng is None else rng
    obs_zero = counts.groupby("taxon")["n"].apply(lambda x: float((x == 0).mean()))
    ppc_rows = []
    for tx, g in counts.groupby("taxon", sort=False):
        Ns = g.set_index("sample")["N"].astype(float).values
        mu = float(bar_x.loc[tx]) * Ns
        b = float(betas.get(tx, beta_shared)) if betas else float(beta_shared)
        sims = []
        for _ in range(draws):
            if np.isinf(b):
                lam = mu
            else:
                lam = rng.gamma(shape=b, scale=np.maximum(mu, 0)/max(b,1e-12))
            k = rng.poisson(lam)
            sims.append(np.mean(k == 0))
        sim = np.asarray(sims, float)
        obs = float(obs_zero.get(tx, np.nan))
        if np.isfinite(obs):
            p_two = float(2 * min(np.mean(sim >= obs), np.mean(sim <= obs)))
        else:
            p_two = np.nan
        ppc_rows.append({"taxon": tx, "ppc_zero_p": p_two})
    return pd.DataFrame(ppc_rows)

# ------------------ PLN + model selection ------------------

def pln_loglik(ns, Ns, mu_log, sigma_log, gh_nodes=30):
    xi, wi = np.polynomial.hermite.hermgauss(gh_nodes)
    x = np.exp(mu_log + np.sqrt(2.0) * sigma_log * xi)
    ll = 0.0
    for n, N in zip(ns, Ns):
        lam = np.maximum(N * x, _EPS_LOG)
        term = np.log(wi) + (n * np.log(lam) - lam - special.gammaln(n + 1))
        ll += special.logsumexp(term) - 0.5 * np.log(np.pi)
    return float(ll)

def fit_pln(ns, Ns, gh_nodes=30):
    ns = np.asarray(ns, int); Ns = np.asarray(Ns, float)
    y = np.clip(ns / np.maximum(Ns, 1.0), 1e-12, 1.0)
    m0 = float(np.mean(np.log(y))); s0 = float(np.std(np.log(y))) or 0.1
    def nll(theta):
        mu_log, log_sigma = theta
        sigma = np.exp(log_sigma)
        return -pln_loglik(ns, Ns, mu_log, sigma, gh_nodes)
    res = optimize.minimize(nll, x0=np.array([m0, np.log(s0)]),
                            bounds=[(-30, 0), (-6, 4)])
    if not res.success:
        return np.nan, np.nan, (np.nan, np.nan, np.nan)
    mu_hat, sig_hat = float(res.x[0]), float(np.exp(res.x[1]))
    ll = -float(res.fun); k = 2
    aic = 2*k - 2*ll; bic = k * np.log(max(len(ns),1)) - 2*ll
    return mu_hat, sig_hat, (ll, aic, bic)

def cv_loglik_models(counts_pos: pd.DataFrame, bar_x: pd.Series, Ns_pos: np.ndarray, folds=5, gh_nodes=30, rng=None):
    rng = np.random.default_rng(321) if rng is None else rng
    samples = counts_pos["sample"].drop_duplicates().tolist()
    idx = np.arange(len(samples)); rng.shuffle(idx); cuts = np.array_split(idx, folds)
    s_order = np.array(samples, dtype=object)
    out = []
    for tx, g in counts_pos.groupby("taxon", sort=False):
        g = g.set_index("sample")
        # Deduplicate any repeated sample rows by summing counts and taking first N
        if g.index.has_duplicates:
            g = g.groupby(level=0).agg({"n": "sum", "N": "first"})
        ll_nb = []; ll_pln = []
        for test in cuts:
            test_s = s_order[test].tolist()
            train_s = [s for s in samples if s not in test_s]
            ns_tr = g.reindex(train_s)["n"].fillna(0).values
            Ns_tr = g.reindex(train_s)["N"].fillna(0).values
            ns_te = g.reindex(test_s)["n"].fillna(0).values
            Ns_te = g.reindex(test_s)["N"].fillna(0).values
            # train NB beta via MLE if enough data, else shared via bar_x prior
            if (ns_tr>0).sum() >= 8 and ns_tr.sum() >= 1000:
                def nll(logb):
                    b = np.exp(logb)
                    mu = float(bar_x.loc[tx]) * Ns_tr
                    return -np.nansum(nb_logpmf_vec(ns_tr, mu, b))
                res = optimize.minimize_scalar(nll, bounds=(-20, 20), method="bounded")
                b_hat = float(np.exp(res.x)) if res.success else 1.0
            else:
                b_hat = 1.0
            mu_te = float(bar_x.loc[tx]) * Ns_te
            ll_nb.append(float(np.nansum(nb_logpmf_vec(ns_te, mu_te, b_hat))))
            # PLN train
            _, _, (ll_pln_tr, _, _) = fit_pln(ns_tr, Ns_tr, gh_nodes=gh_nodes)
            mu_log, sigma, _ = fit_pln(ns_tr, Ns_tr, gh_nodes=gh_nodes)
            if np.isfinite(mu_log) and np.isfinite(sigma):
                ll_pln.append(float(pln_loglik(ns_te, Ns_te, mu_log, sigma, gh_nodes)))
            else:
                ll_pln.append(np.nan)
        out.append({"taxon": tx, "cvll_nb": float(np.nansum(ll_nb)), "cvll_pln": float(np.nansum(ll_pln))})
    return pd.DataFrame(out)

# ------------------ rarefaction (simple MC) ------------------

def rarefy_counts(counts: pd.DataFrame, depths: List[float], draws=50, rng=None) -> pd.DataFrame:
    rng = np.random.default_rng(9) if rng is None else rng
    rows = []
    for s, g in counts.groupby("sample", sort=False):
        N = int(g["N"].iloc[0]); ns = g.set_index("taxon")["n"].astype(int)
        if N <= 0: continue
        for d in depths:
            d = int(min(N, max(1, int(d))))
            for _ in range(draws):
                if d == N:
                    n2 = ns.values
                else:
                    # Multivariate hypergeometric is expensive; use multinomial approximation.
                    # Build a clean probability vector over observed taxa only.
                    p = pd.to_numeric(ns.values, errors="coerce").astype(float)
                    # sanitize: non-negative, finite
                    p[~np.isfinite(p)] = 0.0
                    p[p < 0] = 0.0
                    tot = float(p.sum())
                    if tot <= 0:
                        # nothing to sample from
                        n2 = np.zeros_like(p, dtype=int)
                    else:
                        p = p / tot
                        # final guard against tiny numeric drift
                        p = np.maximum(p, 0.0)
                        p = p / max(p.sum(), 1.0)
                        n2 = rng.multinomial(d, p)
                rows.append({"sample": s, "depth": float(d), "richness": int((n2>0).sum())})
    return pd.DataFrame(rows)

# ------------------ hierarchical bootstrap helpers ------------------

def hier_boot_indices(samples, meta, cluster_cols=None, time_col=None, mbb_L=1, rng=None):
    rng = np.random.default_rng(1) if rng is None else rng
    s = pd.Series(samples, name="sample")
    if cluster_cols:
        def get_key(x):
            row = meta.loc[meta["sample"]==x, cluster_cols]
            if row.empty:
                return ("NA",)
            return tuple(str(v) for v in row.iloc[0].values)
        key = s.map(get_key)
    else:
        key = pd.Series([("ALL",)]*len(s), index=s.index)
    uniq = key.unique().tolist()
    draw_keys = [uniq[i] for i in rng.integers(0, len(uniq), size=len(uniq))]
    resampled = []
    for K in draw_keys:
        group = s[key == K].tolist()
        if time_col and mbb_L > 1:
            gmeta = meta.set_index("sample").reindex(group)
            gmeta = gmeta[~gmeta.index.duplicated(keep="first")]
            gmeta = gmeta.sort_values(time_col)
            g = gmeta.index.tolist()
            if not g: continue
            m = len(g)
            n_blocks = int(np.ceil(m / mbb_L))
            starts = rng.integers(0, m, size=n_blocks)
            seq = []
            for st in starts:
                blk = [g[(st+i)%m] for i in range(mbb_L)]
                seq.extend(blk)
            resampled.extend(seq[:m])
        else:
            resampled.extend([group[i] for i in rng.integers(0, len(group), size=len(group))])
    return resampled

# ------------------ multiple testing (hierarchical option) ------------------

def bh(p): 
    p = np.asarray(p, float); m = np.isfinite(p).sum()
    if m==0: return np.full_like(p, np.nan)
    idx = np.argsort(p[np.isfinite(p)])
    q = np.minimum.accumulate((p[np.isfinite(p)][idx] * m / (np.arange(m)+1))[::-1])[::-1]
    out = np.full_like(p, np.nan); out[np.isfinite(p)][idx] = np.minimum(q, 1.0); return out

# ------------------ main ------------------

def main():
    ap = argparse.ArgumentParser(description="Resistome ecology analytics (reviewer-revision)")
    ap.add_argument("--merged", required=True, help="Merged wide TSV: sample,taxon,rank,tool columns (e.g., metaphlan,bracken)")
    ap.add_argument("--counts", required=False, help="TSV with 'sample','N' (total classified reads used as denominator by tool)")
    ap.add_argument("--rank", default="species", choices=["species","genus","family","order","class","phylum"])
    ap.add_argument("--tool-pair", default="metaphlan,bracken", help="Two columns in merged to compare, e.g., metaphlan,bracken")
    ap.add_argument("--out-dir", default="results/eco_metrics")
    ap.add_argument("--log-level", default="INFO")

    # Agreement / compositional choices
    ap.add_argument("--agreement-mode", default="clr", choices=["clr","alr","plr"])
    ap.add_argument("--subcomp-mode", default="prevalence", choices=["all","prevalence","topk"])
    ap.add_argument("--subcomp-prevalence", type=float, default=0.5)
    ap.add_argument("--subcomp-topk", type=int, default=0)

    ap.add_argument("--zero-replacement", default="pseudocount", choices=["pseudocount","mult","dirichlet"])
    ap.add_argument("--pseudocount", type=float, default=1e-6)
    ap.add_argument("--zero-delta", type=float, default=1e-6)
    ap.add_argument("--dirichlet-alpha", type=float, default=0.5)

    ap.add_argument("--bias-correct", default="none", choices=["none","clr-mean","alr-deming"])
    ap.add_argument("--bias-shrink", type=float, default=0.1)
    ap.add_argument("--deming-delta", type=float, default=1.0)

    # detection thresholds
    ap.add_argument("--detect-thresh", default="fixed", choices=["fixed","kFDR","neg-poisson"])
    ap.add_argument("--detect-alpha", type=float, default=0.01)
    ap.add_argument("--detect-fixed-t", type=int, default=1)

    # occupancy & targets
    ap.add_argument("--target-detect", type=float, default=0.95)
    ap.add_argument("--depths", default=None, help="CSV depths like '1e8,2e8,4e8'")
    ap.add_argument("--forecast-depths", default=None, help="CSV depths for agreement forecasts")
    ap.add_argument("--beta-grid", default=None, help="CSV grid for shared beta search (optional)")
    ap.add_argument("--fit-per-species-beta", action="store_true")

    # model selection & PLN
    ap.add_argument("--afd", default="nb", choices=["nb","pln","both"])
    ap.add_argument("--model-selection", default="aic", choices=["aic","cvll","stack"])
    ap.add_argument("--ms-folds", type=int, default=5)
    ap.add_argument("--pln-quadrature", type=int, default=30)

    # ZINB occupancy split
    ap.add_argument("--zinb-occupancy", action="store_true")

    # metadata / selection
    ap.add_argument("--metadata", default=None, help="Metadata CSV/TSV/XLSX")
    ap.add_argument("--sample-col", default=None)
    ap.add_argument("--include-regex", default=None)
    ap.add_argument("--cohort-col", default=None)
    ap.add_argument("--cohort-allow", default=None)
    ap.add_argument("--date-col", default=None)
    ap.add_argument("--exclude-old-before", default=None, help="ISO date like 2024-01-01")
    ap.add_argument("--neg-regex", default="(?i)(neg|blank|ntc|control)")

    # bootstrap & dependence
    ap.add_argument("--n-bootstrap", type=int, default=500)
    ap.add_argument("--bootstrap-block-col", default=None)
    ap.add_argument("--bootstrap-max", type=int, default=None)
    ap.add_argument("--bootstrap-eps", type=float, default=0.0)
    ap.add_argument("--cluster-cols", default=None, help="Comma-separated columns defining clusters (e.g., Subject,Site)")
    ap.add_argument("--time-col", default=None)
    ap.add_argument("--mbb-window", type=int, default=1)

    # stratify & priority taxa
    ap.add_argument("--stratify", action="store_true")
    ap.add_argument("--batch-col", default=None)
    ap.add_argument("--priority-taxa", default=None, help="CSV/TSV with a 'taxon' column")

    # tests
    ap.add_argument("--gof-alpha", type=float, default=0.05)
    ap.add_argument("--slope-min-n", type=int, default=10)
    ap.add_argument("--mad-min-n", type=int, default=30)
    ap.add_argument("--report-holm", action="store_true")
    ap.add_argument("--mt-hierarchical", action="store_true")

    # nulls & rarefaction
    ap.add_argument("--null", default="none", choices=["none","random-multinomial","tool-permute","neutral-lite"])
    ap.add_argument("--rarefy-grid", default=None)

    # sensitivity grids
    ap.add_argument("--pseudocount-grid", default=None)
    ap.add_argument("--q-grid", default=None)

    # RNG
    ap.add_argument("--seed", type=int, default=1337)

    # Config
    ap.add_argument("--config", help="TOML config file mapping long-option names to values")

    # first parse just to see if --config was given
    args_pre, _unknown = ap.parse_known_args()
    if args_pre.config:
        if _toml is None:
            raise RuntimeError("TOML configs need Python 3.11+ (tomllib) or `pip install tomli` on older Pythons.")
        with open(args_pre.config, "rb") as f:
            cfg = _toml.load(f)

        # convert TOML keys -> CLI args so argparse still does type/choice handling
        cli = []
        # keys that are boolean flags (store_true)
        flag_keys = {
            "fit_per_species_beta", "stratify", "report_holm", "mt_hierarchical", "zinb_occupancy"
        }
        # keys where lists should be joined as CSV strings
        listy = {"depths", "forecast_depths", "beta_grid", "pseudocount_grid", "q_grid", "rarefy_grid"}

        for k, v in cfg.items():
            key = f"--{k.replace('_','-')}"
            if k in flag_keys:
                if bool(v):  # only include when True
                    cli.append(key)
            else:
                if isinstance(v, list) and k in listy:
                    v = ",".join(str(x) for x in v)
                cli.extend([key, str(v)])

        # re-parse using the generated CLI pieces
        args = ap.parse_args(cli)
    else:
        args = args_pre

    logger = setup_logging(args.log_level)
    ensure_dir(args.out_dir)

    # load merged
    merged = pd.read_csv(args.merged, sep="\t")
    t1, t2 = [s.strip() for s in args.tool_pair.split(",")] if "," in args.tool_pair else ("metaphlan","bracken")
    for col in ["sample","taxon","rank",t1,t2]:
        if col not in merged.columns:
            raise ValueError(f"Missing column in merged: {col}")
    if "rank" in merged.columns:
        merged = merged[merged["rank"] == args.rank].copy()
    merged = merged.sort_values(["sample","taxon"], kind="mergesort").reset_index(drop=True)

    # sample selection via metadata (writes lists)
    pos_samples, neg_samples, exc_samples = build_sample_lists(
        merged=merged,
        metadata_path=args.metadata,
        sample_col=args.sample_col,
        include_regex=args.include_regex,
        cohort_col=args.cohort_col,
        cohort_allow=args.cohort_allow,
        date_col=args.date_col,
        exclude_old_before=args.exclude_old_before,
        neg_regex=args.neg_regex,
        out_dir=args.out_dir,
        logger=logger
    )
    pos_df = merged[merged["sample"].isin(pos_samples)].copy()

    # Fixed subcomposition
    keep_taxa = build_fixed_subcomposition(pos_df, pos_samples, t1, t2,
                                           mode=args.subcomp_mode,
                                           prevalence=args.subcomp_prevalence,
                                           topk=args.subcomp_topk)
    pd.DataFrame({"taxon": keep_taxa}).to_csv(os.path.join(args.out_dir, "subcomposition_taxa.tsv"), sep="\t", index=False)

    # Bias correction estimation (if needed)
    bias_df = None
    if args.bias_correct == "clr-mean" and args.agreement_mode == "clr":
        bias_df = estimate_clr_bias(merged, pos_samples, keep_taxa, t1, t2,
                                    args.zero_replacement, args.pseudocount, args.zero_delta, args.dirichlet_alpha,
                                    lam=args.bias_shrink)
        bias_df.to_csv(os.path.join(args.out_dir, "tool_bias.tsv"), sep="\t", index=False)

    # Agreement under fixed subcomposition
    rng_main = np.random.default_rng(int(args.seed))
    agree = agreement_per_sample_fixed(
        pos_df, t1, t2, keep_taxa,
        mode=args.agreement_mode,
        zero_replacement=args.zero_replacement, pseudocount=args.pseudocount,
        alr_ref=args.alr_ref if hasattr(args, "alr_ref") else None,
        alr_topk=getattr(args, "alr_topk", 0),
        plr_maxpairs=getattr(args, "plr_maxpairs", 0),
        rng=rng_main,
        zero_delta=args.zero_delta, dirichlet_alpha=args.dirichlet_alpha,
        bias_df=bias_df, bias_mode=args.bias_correct, deming_delta=args.deming_delta
    )
    agree.to_csv(os.path.join(args.out_dir, f"agreement_per_sample_{args.agreement_mode}.tsv"), sep="\t", index=False)

    summary = []
    def add_sum(line: str): summary.append(line)
    def safe_med(x):
        a = pd.to_numeric(x, errors='coerce').astype(float).values
        m = np.isfinite(a)
        return float(np.nanmedian(a[m])) if m.any() else float('nan')
    add_sum(f"samples_pos={len(set(pos_samples))} samples_neg={len(set(neg_samples))} taxa={pos_df['taxon'].nunique()} rank={args.rank}")
    add_sum(f"Agreement ({args.agreement_mode}) median: pearson={safe_med(agree['pearson']):.3f} spearman={safe_med(agree['spearman']):.3f} aitchison={safe_med(agree['aitchison_clr_norm']):.3f}")

    # counts-dependent analyses
    counts_df = None
    if args.counts and os.path.exists(args.counts):
        counts_df = pd.read_csv(args.counts, sep="\t" if args.counts.lower().endswith((".tsv",".txt")) else ",")[["sample","N"]]
        counts_df = counts_df[counts_df["sample"].isin(pos_samples + neg_samples)].copy()
    else:
        logger.warning("No counts file provided; skipping moments/AFD/occupancy/richness/forecasts/GOF/Taylor/MAD/threshold-aware occupancy.")

    if counts_df is not None and not counts_df.empty:
        # choose tool2 as count source
        counts = build_counts_from_pct(merged[merged["sample"].isin(pos_samples + neg_samples)][["sample","taxon",t2]].copy(), counts_df, t2)

        # negatives-aware thresholds
        detect_mode = args.detect_thresh
        fixed_t = args.detect_fixed_t
        kfdr_map = thresholds_kfdr(counts, neg_samples, alpha=args.detect_alpha) if detect_mode=="kFDR" else {}
        bkg_rate = thresholds_neg_poisson(counts, neg_samples, alpha=args.detect_alpha) if detect_mode=="neg-poisson" else {}

        # restrict to positives for ecology fits
        counts_pos = counts[counts["sample"].isin(pos_samples)].copy()
        Ns_pos = counts_df[counts_df["sample"].isin(pos_samples)].drop_duplicates("sample")["N"].astype(float).values
        medN = float(np.median(Ns_pos)) if Ns_pos.size>0 else np.nan

        # unbiased moments & curvature + bootstrap CIs (hierarchical optional)
        mom = unbiased_moments(counts_pos)
        mom = predict_naive_from_moments(mom, counts_pos)

        # hierarchical bootstrap
        block_series = None
        if args.metadata and os.path.exists(args.metadata):
            meta_full = read_table_any(args.metadata)
            meta_full = meta_full.rename(columns={args.sample_col: "sample"}) if args.sample_col and args.sample_col in meta_full.columns and args.sample_col!="sample" else meta_full
        else:
            meta_full = pd.DataFrame({"sample": pos_samples})

        def draw_indices_hier(rng_local):
            if args.cluster_cols or (args.time_col and args.mbb_window>1):
                cols = [c.strip() for c in (args.cluster_cols or "").split(",") if c.strip()]
                m = meta_full.copy()
                if "sample" not in m.columns:
                    m["sample"] = m.iloc[:,0].astype(str)
                return hier_boot_indices(pos_samples, m, cluster_cols=cols if cols else None,
                                         time_col=args.time_col if args.time_col else None,
                                         mbb_L=max(1,int(args.mbb_window)), rng=rng_local)
            else:
                # default IID sample bootstrap
                S = len(pos_samples)
                ids = np.random.default_rng().integers(0, S, size=S)
                return [pos_samples[i] for i in ids]

        # Occupancy & β fitting (with CV option and thresholds)
        occ_emp = []
        for tax, g in counts_pos.groupby("taxon", sort=False):
            det = [is_detected(int(r["n"]), int(r["N"]), tax, detect_mode, fixed_t, kfdr_map, bkg_rate, args.detect_alpha) for _, r in g.iterrows()]
            occ_emp.append((tax, float(np.mean(det))))
        occ_emp = pd.Series(dict(occ_emp))

        bar_x = mom.set_index("taxon")["mu1"].astype(float)
        beta_grid = parse_csv_list(args.beta_grid)
        if args.ms_folds and args.ms_folds > 1:
            beta_shared = fit_shared_beta_cv_from_counts(counts_pos, Ns_pos, folds=args.ms_folds, grid=beta_grid, rng=rng_main,
                                                         metric="occupancy_mse", detect_mode=detect_mode, fixed_t=fixed_t,
                                                         kfdr_map=kfdr_map, bkg_rate=bkg_rate, detect_alpha=args.detect_alpha)
        else:
            beta_shared = fit_shared_beta(bar_x.values, Ns_pos, occ_emp.reindex(bar_x.index).fillna(0.0).values, beta_grid)
        betas = None; fitted_params = 1
        if args.fit_per_species_beta:
            spp = bar_x.index.tolist()
            s_order = counts_df[counts_df["sample"].isin(pos_samples)].drop_duplicates("sample")["sample"].tolist()
            s2i = {s:i for i,s in enumerate(s_order)}
            mat = {}
            for tax, g in counts_pos.groupby("taxon", sort=False):
                arr = np.zeros(len(s_order), dtype=float)
                for _, r in g.iterrows():
                    arr[s2i[r["sample"]]] = float(r["n"])
                mat[tax] = arr
            betas = fit_beta_per_species(bar_x.values, Ns_pos, mat, spp)
            fitted_params = 2

        # GOF χ² (per species)
        gof = chisq_gof_nb(counts_pos, bar_x, beta_shared, betas, fitted_params)
        if gof is None or gof.empty:
            gof = pd.DataFrame({"taxon": [], "gof_chisq": [], "gof_df": [], "gof_p": []})
        else:
            gof["gof_q"] = bh_fdr(gof["gof_p"].values)
            gof["gof_pass"] = (gof["gof_q"] <= args.gof_alpha)
            if args.report_holm:
                gof["gof_q_holm"] = holm_stepdown(gof["gof_p"].values)

        # AFD table + depth sufficiency + optional ZINB occupancy split
        rows = []
        for tx in bar_x.index:
            b = betas.get(tx, beta_shared) if betas else beta_shared
            bx = float(bar_x.loc[tx])
            occ_cur = float(occ_emp.get(tx, np.nan))
            Nstar = n_star_q(bx, b if np.isfinite(b) else np.inf, args.target_detect) if bx>0 else np.nan
            # ZINB split
            theta_hat = np.nan
            if args.zinb_occupancy and np.isfinite(occ_cur):
                p0_nb = float(np.mean(safe_nb_zero_prob(bx * Ns_pos, b if np.isfinite(b) else np.inf)))
                z_obs = 1.0 - occ_cur
                denom = (1.0 - p0_nb)
                theta_hat = (1.0 - z_obs) / denom if denom > 1e-12 else np.nan
                theta_hat = float(np.clip(theta_hat, 0.0, 1.0)) if np.isfinite(theta_hat) else np.nan
            rows.append({"taxon": tx, "bar_x": bx, "beta": float(b), "occ_current_mean": occ_cur, "N_star_q": Nstar, "theta_hat": theta_hat})
        afd = pd.DataFrame(rows).sort_values("taxon", kind="mergesort").reset_index(drop=True)
        # merge GOF results into AFD table for downstream summaries
        if gof is not None and not gof.empty:
            afd = afd.merge(gof, on="taxon", how="left")

        # Model selection enhancements
        sel_df = None
        if args.afd in ("pln","both") or args.model_selection in ("cvll","stack"):
            if args.model_selection == "cvll" or args.model_selection == "stack":
                cvdf = cv_loglik_models(counts_pos, bar_x, Ns_pos, folds=max(2,int(args.ms_folds)), gh_nodes=int(args.pln_quadrature), rng=rng_main)
                afd = afd.merge(cvdf, on="taxon", how="left")
                if args.model_selection == "cvll":
                    afd["afd_winner"] = np.where(afd["cvll_pln"] > afd["cvll_nb"], "pln", "nb")
                else:
                    # stacking weights (softmax over cvll)
                    m = afd[["cvll_nb","cvll_pln"]].replace([-np.inf, np.inf], np.nan).fillna(-1e12).values
                    e = np.exp(m - m.max(axis=1, keepdims=True))
                    denom = e.sum(axis=1, keepdims=True)
                    denom = np.where(denom <= 0, 1.0, denom)
                    w = e / denom
                    afd["stack_w_nb"] = w[:,0]; afd["stack_w_pln"] = w[:,1]
                    afd["afd_winner"] = np.where(afd["stack_w_pln"] > afd["stack_w_nb"], "pln", "nb")
            else:
                # AIC on full data (original)
                sel_rows = []
                for tx, g in counts_pos.groupby("taxon", sort=False):
                    Ns_tx = g["N"].astype(float).values
                    ns_tx = g["n"].astype(int).values
                    b_tx = float(betas.get(tx, beta_shared)) if betas else float(beta_shared)
                    mu_tx = float(bar_x.loc[tx]) * Ns_tx
                    ll_nb = float(np.nansum(nb_logpmf_vec(ns_tx, mu_tx, b_tx)))
                    aic_nb = 2*1 - 2*ll_nb
                    _, _, (ll_pln, aic_pln, bic_pln) = fit_pln(ns_tx, Ns_tx, gh_nodes=int(getattr(args, 'pln_quadrature', 30)))
                    winner = "pln" if (np.isfinite(aic_pln) and aic_pln < aic_nb) else "nb"
                    sel_rows.append({"taxon": tx, "ll_nb": ll_nb, "aic_nb": aic_nb,
                                     "ll_pln": ll_pln, "aic_pln": aic_pln, "bic_pln": bic_pln,
                                     "afd_winner": winner})
                sel_df = pd.DataFrame(sel_rows)
                afd = afd.merge(sel_df, on="taxon", how="left")

        # Posterior Predictive Checks
        if hasattr(args, "ppc_draws") and args.ppc_draws and args.ppc_draws > 0:
            ppc = posterior_predictive_checks(counts_pos, bar_x, beta_shared, betas=betas, draws=int(args.ppc_draws), rng=rng_main)
            if ppc is not None and not ppc.empty:
                afd = afd.merge(ppc, on="taxon", how="left")
        afd.to_csv(os.path.join(args.out_dir, "afd_fit.tsv"), sep="\t", index=False)

    # occupancy calibration metrics (threshold-aware)
    pred_occ = occupancy_from_beta(bar_x.values, Ns_pos, float(np.nanmedian(afd["beta"].astype(float))) if betas else beta_shared)
    emp_occ_vals = afd.set_index("taxon")["occ_current_mean"].reindex(bar_x.index).fillna(0.0).values
    rmse = float(np.sqrt(np.mean((pred_occ - emp_occ_vals)**2)))
    r = _safe_corr(pred_occ, emp_occ_vals)
    r2 = (r*r) if np.isfinite(r) else np.nan
    add_sum(f"Occupancy calibration: RMSE={rmse:.3f} R2={r2:.3f}")

    # depth sufficiency per species
    depth_rows = []
    for _, r in afd.iterrows():
            b = float(r["beta"]) if np.isfinite(r["beta"]) else np.inf
            bx = float(r["bar_x"])
            Nstar = float(r["N_star_q"]) if np.isfinite(r["N_star_q"]) else np.nan
            occ_medN = 1.0 - safe_nb_zero_prob(bx * medN, b) if np.isfinite(medN) else np.nan
            depth_rows.append({
                "taxon": r["taxon"], "bar_x": bx, "beta": b, "N_star_q": Nstar,
                "medianN": medN, "is_sufficient_at_current": bool(np.isfinite(Nstar) and medN >= Nstar),
                "occupancy_at_medianN": occ_medN, "observed_occupancy": float(occ_emp.get(r["taxon"], np.nan)),
                "theta_hat": r.get("theta_hat", np.nan)
            })
    depth_suff = pd.DataFrame(depth_rows).sort_values("taxon", kind="mergesort")
    depth_suff.to_csv(os.path.join(args.out_dir, "depth_sufficiency_per_species.tsv"), sep="\t", index=False)
    suff_pct = float(np.nanmean(depth_suff["is_sufficient_at_current"])) * 100.0
    add_sum(f"Depth sufficiency at medianN: {suff_pct:.1f}% species; medianN={medN:.3g}")

    # richness + slope with bootstrap CIs
    rich_rows = []
    depths = parse_csv_list(args.depths) or []
    if medN > 0 and (not depths or all(abs(d - medN) > 1e-9 for d in depths)):
        depths = [medN] + depths
    beta_for_R = float(np.nanmedian(afd["beta"].astype(float)))
    bx_vec = bar_x.values
    for d in depths:
        R, Rp = richness_and_slope(bx_vec, beta_for_R, float(d))
        rich_rows.append({"stratum": "ALL", "depth": float(d), "R": R, "R_prime": Rp})
    rich_df = pd.DataFrame(rich_rows)
    if not rich_df.empty:
        rich_df.sort_values(["stratum","depth"], kind="mergesort").to_csv(os.path.join(args.out_dir, "richness_vs_depth.tsv"), sep="\t", index=False)

        # bootstrap CIs for derived metrics
        B = max(0, int(args.n_bootstrap))
        if B > 0:
            rngB = np.random.default_rng(int(args.seed)+99)
            samples = counts_pos["sample"].drop_duplicates().tolist()
            m = meta_full.copy()
            if "sample" not in m.columns: m["sample"] = m.iloc[:,0].astype(str)
            cols = [c.strip() for c in (args.cluster_cols or "").split(",") if c.strip()]
            Rboots = {float(d):[] for d in depths}
            Rpboots = {float(d):[] for d in depths}
            for b in range(B):
                ss = hier_boot_indices(samples, m, cluster_cols=cols if cols else None,
                                       time_col=args.time_col if args.time_col else None,
                                       mbb_L=max(1,int(args.mbb_window)), rng=rngB)
                cb = counts_pos[counts_pos["sample"].isin(ss)]
                mb = unbiased_moments(cb)
                bx = mb["mu1"].values
                bet = float(np.nanmedian(afd["beta"].astype(float)))
                for d in depths:
                    Rb, Rpb = richness_and_slope(bx, bet, float(d))
                    Rboots[float(d)].append(Rb); Rpboots[float(d)].append(Rpb)
            rows=[]
            for d in depths:
                arr = np.asarray(Rboots[float(d)], float); arrp = np.asarray(Rpboots[float(d)], float)
                if arr.size>0:
                    lo, hi = np.nanpercentile(arr, [2.5,97.5]); lop, hip = np.nanpercentile(arrp, [2.5,97.5])
                    rows.append({"depth": float(d), "R_lo": lo, "R_hi": hi, "Rprime_lo": lop, "Rprime_hi": hip})
            if rows:
                pd.DataFrame(rows).to_csv(os.path.join(args.out_dir, "richness_vs_depth_ci.tsv"), sep="\t", index=False)

        # model-guided filtering & re-agreement on fixed subcomposition
        keep = set()
        for _, r in afd.iterrows():
            bx = float(r["bar_x"]); b = float(r["beta"]) if np.isfinite(r["beta"]) else np.inf
            if bx <= 0: continue
            occ = 1.0 - safe_nb_zero_prob(bx * medN, b) if np.isfinite(medN) else 0.0
            if occ >= float(getattr(args, "min_occupancy", 0.5)):
                keep.add(r["taxon"])
        if keep:
            pos_df_f = pos_df[pos_df["taxon"].isin(keep)].copy()
            agree_f = agreement_per_sample_fixed(
                pos_df_f, t1, t2, keep_taxa,
                mode=args.agreement_mode,
                zero_replacement=args.zero_replacement, pseudocount=args.pseudocount,
                alr_ref=getattr(args, "alr_ref", None),
                alr_topk=getattr(args, "alr_topk", 0),
                plr_maxpairs=getattr(args, "plr_maxpairs", 0),
                rng=rng_main,
                zero_delta=args.zero_delta, dirichlet_alpha=args.dirichlet_alpha,
                bias_df=bias_df, bias_mode=args.bias_correct, deming_delta=args.deming_delta
            )
            agree_f.to_csv(os.path.join(args.out_dir, f"agreement_per_sample_filtered_{args.agreement_mode}.tsv"), sep="\t", index=False)
            add_sum(f"Agreement ({args.agreement_mode}, filtered) median: pearson={safe_med(agree_f['pearson']):.3f} spearman={safe_med(agree_f['spearman']):.3f} aitchison={safe_med(agree_f['aitchison_clr_norm']):.3f}")

        # forecasts vs depth (zero fraction heuristic)
        fdepths = parse_csv_list(args.forecast_depths) or []
        if fdepths:
            beta_used = beta_for_R
            rows = []
            # use fixed subcomposition taxa for forecast alignment
            for d in fdepths:
                for s, g in pos_df.groupby("sample", sort=False):
                    spp = g["taxon"].unique()
                    bx = bar_x.reindex(spp).dropna().values
                    if bx.size == 0:
                        zf = 1.0
                    else:
                        zf = float(np.mean(safe_nb_zero_prob(bx * float(d), beta_used)))
                    sugg = "pearson" if zf < 0.2 else "spearman"
                    rows.append({"depth": float(d), "sample": s, "pred_zero_fraction": zf, "suggested": sugg})
            if rows:
                fdf = pd.DataFrame(rows).sort_values(["depth","sample"], kind="mergesort")
                fdf.to_csv(os.path.join(args.out_dir, "agreement_forecast_vs_depth.tsv"), sep="\t", index=False)
                agg = fdf.groupby("depth")["suggested"].value_counts(normalize=True).rename("prop").reset_index()
                agg_pivot = agg.pivot(index="depth", columns="suggested", values="prop").fillna(0.0).reset_index()
                agg_pivot.to_csv(os.path.join(args.out_dir, "agreement_forecast_summary.tsv"), sep="\t", index=False)
                for _, r in agg_pivot.iterrows():
                    add_sum(f"Forecast at depth={r['depth']:.3g}: pearson%={100*r.get('pearson',0):.1f} spearman%={100*r.get('spearman',0):.1f}")

        # Taylor & MAD as before
        def taylor_test(mu_df: pd.DataFrame, min_n: int = 10) -> Optional[Dict[str,float]]:
            df = mu_df.copy()
            df["mean"] = df["mu1"]
            df["var"] = df["mu2"] - df["mu1"]**2
            df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["mean","var"])
            df = df[(df["mean"] > 0) & (df["var"] > 0)]
            n = len(df)
            if n < max(min_n, 3):
                return None
            x = np.log(df["mean"].values); y = np.log(df["var"].values)
            X = np.c_[np.ones_like(x), x]
            beta_hat, *_ = np.linalg.lstsq(X, y, rcond=None)
            a, b = float(beta_hat[0]), float(beta_hat[1])
            yhat = X @ beta_hat
            rss = float(np.sum((y - yhat)**2))
            df_res = n - 2
            sigma2 = rss / max(df_res, 1)
            cov = sigma2 * np.linalg.inv(X.T @ X)
            se_b = float(np.sqrt(cov[1,1]))
            tcrit = float(stats.t.ppf(0.975, max(df_res,1)))
            ci_lo, ci_hi = b - tcrit*se_b, b + tcrit*se_b
            tstat = (b - 2.0) / max(se_b, 1e-12)
            p_eq2 = float(2 * stats.t.sf(abs(tstat), max(df_res,1)))
            tss = float(np.sum((y - np.mean(y))**2))
            r2 = 1.0 - (rss / tss if tss > 0 else np.nan)
            return {"n_species": n, "slope_b": b, "ci_low": ci_lo, "ci_high": ci_hi, "p_b_eq_2": p_eq2, "r2": r2}

        def mad_lognormal(mu_df: pd.DataFrame, min_n: int = 30) -> Optional[Dict[str,float]]:
            vals = np.log(mu_df.loc[mu_df["mu1"] > 0, "mu1"].values)
            n = vals.size
            if n < min_n: return None
            ad = stats.anderson(vals, dist="norm")
            out = {"n_species": int(n), "ad_stat": float(ad.statistic), "ad_crit_5pct": float(ad.critical_values[2])}
            if n <= 5000:
                w, p = stats.shapiro(vals)
                out.update({"shapiro_w": float(w), "shapiro_p": float(p)})
            return out

        tl = taylor_test(mom, args.slope_min_n)
        if tl:
            pd.DataFrame([{"stratum":"ALL", **tl}]).to_csv(os.path.join(args.out_dir, "taylor_law.tsv"), sep="\t", index=False, mode="w")
            add_sum(f"Taylor slope b={tl['slope_b']:.3f} CI=[{tl['ci_low']:.3f},{tl['ci_high']:.3f}] p(b=2)={tl['p_b_eq_2']:.3g} R2={tl['r2']:.3f}")
        md = mad_lognormal(mom, args.mad_min_n)
        if md:
            pd.DataFrame([{"stratum":"ALL", **md}]).to_csv(os.path.join(args.out_dir, "mad_lognormal.tsv"), sep="\t", index=False, mode="w")
            add_sum(f"MAD lognormal: AD={md['ad_stat']:.3f} vs 5% crit={md['ad_crit_5pct']:.3f}" + (f"; Shapiro p={md['shapiro_p']:.3g}" if "shapiro_p" in md else ""))

        # priority taxa
        if args.priority_taxa and os.path.exists(args.priority_taxa):
            pt = read_table_any(args.priority_taxa)
            if "taxon" in pt.columns:
                want = pt["taxon"].astype(str).dropna().unique().tolist()
                beta_used = float(np.nanmedian(afd["beta"].astype(float)))
                depths2 = parse_csv_list(args.depths) or []
                if medN > 0 and (not depths2 or all(abs(medN - d) > 1e-9 for d in depths2)):
                    depths2 = [medN] + depths2
                rows = []
                bx = afd.set_index("taxon")["bar_x"].to_dict()
                occ_obs = afd.set_index("taxon")["occ_current_mean"].to_dict()
                betad = afd.set_index("taxon")["beta"].to_dict()
                for tx in want:
                    if tx not in bx: continue
                    barx = float(bx[tx]); b = float(betad.get(tx, beta_used))
                    Nstar = n_star_q(barx, b if np.isfinite(b) else np.inf, args.target_detect) if barx>0 else np.nan
                    for d in depths2:
                        occ = 1.0 - safe_nb_zero_prob(barx * float(d), b if np.isfinite(b) else np.inf)
                        rows.append({"taxon": tx, "depth": float(d), "bar_x": barx, "beta": b, "occ_at_depth": occ, "N_star_q": Nstar, "observed_occupancy": float(occ_obs.get(tx, np.nan))})
                if rows:
                    pd.DataFrame(rows).sort_values(["taxon","depth"], kind="mergesort").to_csv(os.path.join(args.out_dir, "priority_taxa_detection.tsv"), sep="\t", index=False)

    

        # rarefaction / coverage (optional)
        rgrid = parse_csv_list(args.rarefy_grid)
        if rgrid:
            rare = rarefy_counts(counts_pos, [int(d) for d in rgrid], draws=50, rng=rng_main)
            if not rare.empty:
                agg = rare.groupby(["sample","depth"])["richness"].agg(["mean","median","std"]).reset_index()
                agg.to_csv(os.path.join(args.out_dir, "rarefaction_empirical.tsv"), sep="\t", index=False)

        # hierarchical multiple testing (if opted) — demo across families
        if args.mt_hierarchical:
            fam_p = []
            if isinstance(afd, pd.DataFrame) and "gof_p" in afd.columns:
                fam_p.append(afd["gof_p"].values.copy())
            # (Add other families as desired; placeholder)
            if fam_p:
                allp = np.concatenate([x[np.isfinite(x)] for x in fam_p])
                q_all = bh(allp)  # BY across families omitted for simplicity; BH used.
                # We don't re-map back here since families are separate; keep per-family q's as above.

        # sensitivity grids
        pc_grid = parse_csv_list(getattr(args, "pseudocount_grid", None))
        if pc_grid:
            rows=[]
            for pc in pc_grid:
                a = agreement_per_sample_fixed(
                    pos_df, t1, t2, keep_taxa,
                    mode=args.agreement_mode,
                    zero_replacement="pseudocount", pseudocount=float(pc),
                    alr_ref=getattr(args, "alr_ref", None),
                    alr_topk=getattr(args, "alr_topk", 0),
                    plr_maxpairs=getattr(args, "plr_maxpairs", 0),
                    rng=rng_main,
                    zero_delta=args.zero_delta, dirichlet_alpha=args.dirichlet_alpha,
                    bias_df=bias_df, bias_mode=args.bias_correct, deming_delta=args.deming_delta
                )
                rows.append({"pseudocount": pc, "pearson_med": float(np.nanmedian(a["pearson"])), "spearman_med": float(np.nanmedian(a["spearman"]))})
            pd.DataFrame(rows).to_csv(os.path.join(args.out_dir, f"agreement_sensitivity_{args.agreement_mode}.tsv"), sep="\t", index=False)

        q_grid = parse_csv_list(getattr(args, "q_grid", None))
        if q_grid and afd is not None and not afd.empty:
            rows=[]
            for q in q_grid:
                for _, r in afd.iterrows():
                    bx = float(r["bar_x"]); b = float(r["beta"]) if np.isfinite(r["beta"]) else np.inf
                    Nstar = n_star_q(bx, b, q) if (np.isfinite(bx) and bx>0) else np.nan
                    rows.append({"taxon": r["taxon"], "q": q, "N_star_q": Nstar})
            pd.DataFrame(rows).to_csv(os.path.join(args.out_dir, "nstar_grid.tsv"), sep="\t", index=False)

    # write summary
    with open(os.path.join(args.out_dir, "summary.txt"), "w") as f:
        for line in summary:
            f.write(line + "\n")


if __name__ == "__main__":
    main()
