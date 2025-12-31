#!/usr/bin/env python3
# Portfolio-safe copy: compute alpha/beta diversity from Bracken outputs (species-level).
# Identifiers generalized; no private data included.

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from skbio.diversity import alpha_diversity, beta_diversity
except Exception as e:
    raise SystemExit("This script requires scikit-bio. Install with: pip install scikit-bio") from e


def read_bracken_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    # try common bracken cols
    if "name" not in df.columns:
        raise ValueError(f"Bracken file missing 'name' column: {path}")
    if "new_est_reads" in df.columns:
        valcol = "new_est_reads"
    elif "est_reads" in df.columns:
        valcol = "est_reads"
    else:
        raise ValueError(f"Bracken file missing abundance column: {path}")
    return df[["name", valcol]].rename(columns={valcol: "count"})


def build_count_matrix(bracken_files: list[str]) -> tuple[pd.DataFrame, list[str]]:
    tables = []
    sample_ids = []
    for fp in bracken_files:
        sid = os.path.basename(fp).replace(".bracken", "")
        sample_ids.append(sid)
        t = read_bracken_table(fp).set_index("name")
        t.columns = [sid]
        tables.append(t)
    mat = pd.concat(tables, axis=1).fillna(0.0)
    # enforce numeric
    mat = mat.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return mat, sample_ids


def compute_alpha(mat: pd.DataFrame, metrics: list[str]) -> pd.DataFrame:
    counts = mat.T.values  # samples x features
    ids = mat.columns.tolist()
    out = {}
    for m in metrics:
        vals = alpha_diversity(m, counts, ids=mat.T.index.astype(str))
        out[m] = vals
    return pd.DataFrame(out)


def compute_beta(mat: pd.DataFrame, metrics: list[str]) -> dict[str, "skbio.stats.distance.DistanceMatrix"]:
    counts = mat.T.values
    ids = mat.columns.tolist()
    dm = {}
    for m in metrics:
        dm[m] = beta_diversity(m, counts, ids=mat.T.index.astype(str))
    return dm


def main():
    ap = argparse.ArgumentParser(description="Alpha/beta diversity from Bracken outputs.")
    ap.add_argument("bracken_dir", help="Directory containing *.bracken")
    ap.add_argument("--alpha", nargs="+", default=["shannon", "simpson"], help="Alpha metrics for scikit-bio")
    ap.add_argument("--beta", nargs="+", default=["braycurtis", "jaccard"], help="Beta metrics for scikit-bio")
    ap.add_argument("--outdir", default="results/diversity", help="Output directory")
    args = ap.parse_args()

    bracken_files = sorted(glob.glob(os.path.join(args.bracken_dir, "*.bracken")))
    if not bracken_files:
        raise SystemExit(f"No .bracken files found in {args.bracken_dir}")

    mat, sample_ids = build_count_matrix(bracken_files)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Save count matrix
    mat.to_csv(outdir / "counts_wide.tsv", sep="\t")

    # Alpha
    alpha_df = compute_alpha(mat, args.alpha)
    alpha_df.to_csv(outdir / "alpha_diversity.tsv", sep="\t", index=True)

    # Beta
    beta = compute_beta(mat, args.beta)
    for m, dm in beta.items():
        dm.to_data_frame().to_csv(outdir / f"beta_{m}.tsv", sep="\t", index=True)

    print(f"Wrote outputs to: {outdir}")


if __name__ == "__main__":
    main()
