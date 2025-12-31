#!/usr/bin/env python3
"""io_utils.py (portfolio scaffold)

This file is *new* scaffolding to make the comparison modules runnable in a public portfolio.
It replaces private/internal helpers with small, dependency-light equivalents.

Expected input:
- A "joined" table with at least:
  sampleID, species, bracken_fraction, rel_abundance, new_est_reads (optional in some modules)

By default, `read_joined_table()` reads:
  results/tables/joined_metaphlan_bracken.tsv

You can override with env var:
  META_ECO_JOINED_TABLE=/path/to/file.tsv
"""
from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple

DEFAULT_JOINED = "results/tables/joined_metaphlan_bracken.tsv"


def float_or_nan(x: Any) -> float:
    if x is None:
        return math.nan
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x).strip()
    if s == "" or s.lower() in {"na", "nan", "none", "null"}:
        return math.nan
    try:
        return float(s)
    except ValueError:
        return math.nan


def read_tsv(path: str) -> List[Dict[str, str]]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(str(p))
    with p.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [dict(r) for r in reader]


def write_tsv(path: str, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in rows:
            out = {}
            for f in fieldnames:
                v = r.get(f, "")
                if isinstance(v, float) and math.isnan(v):
                    out[f] = ""
                else:
                    out[f] = v
            writer.writerow(out)


def group_by(rows: Iterable[Dict[str, Any]], key: str) -> Dict[str, List[Dict[str, Any]]]:
    out: Dict[str, List[Dict[str, Any]]] = {}
    for r in rows:
        out.setdefault(str(r.get(key, "")), []).append(r)
    return out


def read_joined_table(path: str | None = None) -> List[Dict[str, Any]]:
    joined_path = path or os.environ.get("META_ECO_JOINED_TABLE") or DEFAULT_JOINED
    rows = read_tsv(joined_path)
    # normalize key columns
    out = []
    for r in rows:
        rr = dict(r)
        # coerce common numeric fields if present
        for col in ("bracken_fraction", "rel_abundance", "new_est_reads", "est_reads"):
            if col in rr:
                rr[col] = float_or_nan(rr[col])
        out.append(rr)
    return out


def spearman(x: Sequence[float], y: Sequence[float]) -> float:
    """Spearman rho without scipy (ties handled by average ranks)."""
    if len(x) != len(y) or len(x) < 3:
        return math.nan
    pairs = [(a, b) for a, b in zip(x, y) if not (math.isnan(a) or math.isnan(b))]
    if len(pairs) < 3:
        return math.nan
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]

    def rankdata(v: List[float]) -> List[float]:
        idx = sorted(range(len(v)), key=lambda i: v[i])
        ranks = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[idx[j + 1]] == v[idx[i]]:
                j += 1
            avg_rank = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                ranks[idx[k]] = avg_rank
            i = j + 1
        return ranks

    rx = rankdata(xs)
    ry = rankdata(ys)
    mx = sum(rx) / len(rx)
    my = sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    denx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    deny = math.sqrt(sum((b - my) ** 2 for b in ry))
    if denx == 0 or deny == 0:
        return math.nan
    return num / (denx * deny)


def clr_transform(vals: Sequence[float], pseudocount: float = 1e-12) -> List[float]:
    """Centered log-ratio transform for positive values."""
    v = [float(x) for x in vals if x is not None and not math.isnan(float(x)) and float(x) > 0]
    if not v:
        return []
    v = [x + pseudocount for x in v]
    gm = math.exp(sum(math.log(x) for x in v) / len(v))
    return [math.log(x / gm) for x in v]
