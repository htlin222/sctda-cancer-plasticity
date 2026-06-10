# Cautionary/methods paper — reproducibility

A control framework for topological data analysis of perturbation scRNA-seq, with an
EGFR-mutant lung-cancer case study. All analyses run in the conda `sctda` environment
(`ripser` 0.6, `dreimac` 0.3, `scanpy` 1.11) from raw GEO data already staged under `data/`.

## Files
- `manuscript.md` — full manuscript (title, abstract, contributions, results, methods, discussion)
- `cover_letter.md` — cover letter (target: Bioinformatics)
- figures in `../../figures/cautionary/`

## Analysis → output map

| Script | Produces | Used in |
|--------|----------|---------|
| `scripts/16_lineage_circular_coords.py` | circular coordinates, Rayleigh R, clonal angular spread | Control 1; Fig 2a |
| `scripts/17_dispersion_null_and_stability.py` | covariance-matched Gaussian null, bootstrap CIs | Control 2; Fig 1 |
| `scripts/18_structure_cell_biology.py` | leave-one-cluster-out + DE of structure-driving cells | Control 3 |
| `scripts/19_is_structure_cellcycle.py` | per-timepoint cell-cycle regression check | Control 3 |
| `scripts/21_pooled_cc_decisive.py` | pooled embedding, std vs S/G2M-regressed, bootstrap | Control 3; Fig 2b |
| `scripts/22_clonal_memory_decay.py` | clonal memory M(t) + robustness grid | clonal statistic; Fig 3 |
| `scripts/23_memory_depth_matched_null.py` | depth-stratified null for M(t) | clonal statistic |
| `scripts/25_synthetic_benchmark.py` | synthetic ground-truth benchmark (+ 2nd null) | sensitivity/specificity; Fig 4 |
| `scripts/24_cautionary_figures.py` | Fig 1–4 PDFs | all figures |

PDX cell-cycle robustness (residual excess survives removing cycling cells; G1-only ratio
2.07 [1.37,2.61], 100% of bootstraps) is computed inline from the PDX `.raw` layer.

## Regenerate figures
```bash
conda activate sctda
python scripts/17_dispersion_null_and_stability.py   # -> results/foundation/dispersion_null_summary.json
python scripts/16_lineage_circular_coords.py
python scripts/25_synthetic_benchmark.py             # -> results/foundation/synthetic_benchmark.json
python scripts/24_cautionary_figures.py              # -> figures/cautionary/fig1-4.pdf
```

## Key result JSONs (`results/foundation/`)
`dispersion_null_summary.json`, `circular_coord_summary.json`, `cellcycle_structure_test.csv`,
`clonal_memory_decay.json`, `memory_depth_matched.json`, `synthetic_benchmark.json`.
