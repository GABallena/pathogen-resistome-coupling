#!/usr/bin/env python3
# Portfolio-safe copy: paths/identifiers generalized; inputs not included.

import argparse, csv, glob, os, re, sys
from typing import Dict, List

def norm(s: str) -> str:
    return re.sub(r'[^a-z0-9]+', '_', s.strip().lower())

CAND = {
    "contig": ["contig", "contig_id", "contig name", "sequence name", "orf_from"],
    "orf_id": ["orf_id", "orf id", "orf", "orf_identifier"],
    "aro": ["aro", "aro accession", "aro_name", "aro_name_accession"],
    "arg_class": ["drug_class", "drug class"],
    "arg_family": ["amr_gene_family", "amr gene family", "resistance mechanism", "resistance gene family"],
    "model_type": ["model_type", "model type"],
    "perc_identity": ["%_identity", "perc_identity", "percentage_identity", "pct_identity", "identity"],
    "orf_length_aa": ["orf_length_aa", "orf length (aa)", "orf length", "orf_length"]
}

def pick_indices(header: List[str]) -> Dict[str, int]:
    hmap = {norm(h): i for i, h in enumerate(header)}
    out = {}
    for key, opts in CAND.items():
        for cand in opts:
            if cand in hmap:
                out[key] = hmap[cand]
                break
    return out

def main():
    ap = argparse.ArgumentParser(description="Parse RGI TXT folder into tidy TSV.")
    ap.add_argument("--rgi_dir", default="rgi_output", help="Folder with *_rgi.txt")
    ap.add_argument("--out", default="results/tables/rgi_orfs.tsv")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    files = sorted(glob.glob(os.path.join(args.rgi_dir, "*_rgi.txt")))
    if not files:
        print("No RGI TXT files found in rgi_output/", file=sys.stderr)

    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["sample","contig_id","orf_id","aro","arg_class","arg_family","model_type","perc_identity","orf_length_aa"])
        kept = 0
        for tf in files:
            sample = os.path.basename(tf).replace("_rgi.txt","")
            with open(tf, newline="") as fh:
                header = fh.readline().rstrip("\n").split("\t")
                if not header or len(header) < 2:
                    continue
                idx = pick_indices(header)
                for line in fh:
                    row = line.rstrip("\n").split("\t")
                    if not row or len(row) < 2: 
                        continue
                    contig = row[idx["contig"]] if "contig" in idx and idx["contig"] < len(row) else ""
                    orf_id = row[idx["orf_id"]] if "orf_id" in idx and idx["orf_id"] < len(row) else ""
                    aro = row[idx["aro"]] if "aro" in idx and idx["aro"] < len(row) else ""
                    arg_class = row[idx["arg_class"]] if "arg_class" in idx and idx["arg_class"] < len(row) else ""
                    arg_family = row[idx["arg_family"]] if "arg_family" in idx and idx["arg_family"] < len(row) else ""
                    model_type = row[idx["model_type"]] if "model_type" in idx and idx["model_type"] < len(row) else ""
                    perc_identity = row[idx["perc_identity"]] if "perc_identity" in idx and idx["perc_identity"] < len(row) else ""
                    orf_len = row[idx["orf_length_aa"]] if "orf_length_aa" in idx and idx["orf_length_aa"] < len(row) else ""
                    w.writerow([sample, contig, orf_id, aro, arg_class, arg_family, model_type, perc_identity, orf_len])
                    kept += 1
    print(f"Wrote {kept} rows to {args.out}")

if __name__ == "__main__":
    main()
