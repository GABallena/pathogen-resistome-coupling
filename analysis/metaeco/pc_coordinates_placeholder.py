#!/usr/bin/env python3
# Portfolio-safe copy: identifiers generalized; runnable with included io_utils scaffold.

"""Placeholder generating fake PCA coordinates until real ordination implemented.
Consumes the wide CLR matrix for each tool and emits first two PCs as zeros.
Outputs:
 - comparisons/pca_clr_bracken.tsv
 - comparisons/pca_clr_metaphlan.tsv
"""
from __future__ import annotations
from .io_utils import read_tsv, write_tsv

IN_BR = "results/comparisons/ordination_clr_bracken.tsv"
IN_MP = "results/comparisons/ordination_clr_metaphlan.tsv"
OUT_BR = "results/comparisons/pca_clr_bracken.tsv"
OUT_MP = "results/comparisons/pca_clr_metaphlan.tsv"

def main():
    try:
        br = read_tsv(IN_BR)
        mp = read_tsv(IN_MP)
    except FileNotFoundError:
        return
    out_br=[{"sampleID": r['sampleID'], "PC1": 0.0, "PC2": 0.0} for r in br]
    out_mp=[{"sampleID": r['sampleID'], "PC1": 0.0, "PC2": 0.0} for r in mp]
    write_tsv(OUT_BR, out_br, ["sampleID","PC1","PC2"])
    write_tsv(OUT_MP, out_mp, ["sampleID","PC1","PC2"])

if __name__ == '__main__':
    main()
