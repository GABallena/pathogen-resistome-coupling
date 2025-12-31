# pathogen-resistome-coupling

Analysis templates for quantifying and visualizing **pathogen–resistome relationships** from metagenomics-derived tables.

The goal is to make “coupling” analyses reproducible:
- align taxa abundance with AMR feature abundance
- compute correlations / associations robustly
- visualize patterns with clear assumptions and QC checks

## Typical workflow
1. Import/standardize:
   - taxa tables (e.g., genus/species abundance)
   - resistome tables (e.g., ARG family/class abundance)
   - sample metadata (site/time/condition)
2. QC + filtering:
   - library size checks, prevalence thresholds
   - compositional considerations (transformations)
3. Association analyses:
   - correlation / regression / permutation-based checks
4. Visualization:
   - heatmaps, network-style summaries, ranked associations

## Inputs
This repo assumes you provide your own `data/` tables locally.

Recommended structure:
```
data/
  metadata.csv
  taxa_table.tsv
  resistome_table.tsv
```

## License
MIT — see `LICENSE`.
