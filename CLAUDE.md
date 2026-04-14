# scTDA Cancer Plasticity

> Topological data analysis reveals cyclic cell-state plasticity underlying drug tolerance in EGFR-mutant lung cancer.

## What this is

A publication-grade computational biology project applying persistent homology and Mapper algorithms to single-cell RNA-seq data from drug-treated EGFR-mutant NSCLC models. The goal is to discover non-tree-like (cyclic) topological structures in cell-state spaces that standard trajectory inference methods miss, and to link these structures to drug tolerance mechanisms. Target journals: Nature Communications, Genome Biology, Cell Reports.

## Stack

- **Language:** Python 3.11+
- **Environment:** conda (mamba preferred)
- **Core scientific:** scanpy, anndata, ripser, persim, giotto-tda, kepler-mapper
- **Standard scRNA-seq:** scanpy, scvi-tools, leidenalg, bbknn
- **Visualization:** matplotlib, seaborn, plotly
- **Statistics:** scipy, statsmodels, scikit-learn
- **Manuscript figures:** matplotlib with publication-quality settings (300 dpi, vector PDF)
- **Test runner:** pytest
- **Linter:** ruff

## Key commands

```bash
# Create environment
mamba env create -f envs/environment.yml
conda activate sctda

# Run full pipeline
make pipeline

# Run specific phases
make download    # Phase 1: download GEO data
make preprocess  # Phase 2: QC, normalize, HVG, PCA
make topology    # Phase 3: persistent homology + Mapper
make biology     # Phase 4: gene connectivity + pathway enrichment
make validate    # Phase 5: validation on independent dataset
make figures     # Phase 6: publication figures

# Run tests
make test

# Lint
make lint
```

## Project conventions

- Notebook naming: `NN-description.ipynb` (e.g., `01-eda.ipynb`)
- Figure naming: `fig{N}_{panel}_{description}.pdf` (e.g., `fig2_a_persistence_barcode.pdf`)
- Supplementary: `sfig{N}_{description}.pdf`
- Data files: never committed to git; downloaded via `make download`
- Branch naming: `feat/`, `fix/`, `analysis/`
- Commit style: Conventional Commits (`feat:`, `fix:`, `analysis:`, `fig:`)

## Important constraints

- **No raw data in repo** — all data downloaded from GEO/Zenodo via scripts
- **Reproducibility** — every figure must be regenerable from `make figures`
- **Random seeds** — all stochastic operations seeded (seed=42 default)
- **Statistical rigor** — every topological claim requires permutation test (n≥500)
- **Cell cycle ablation** — every H₁ finding must be tested with/without cell cycle genes

## Current focus

Phase 1 — see `plan/01-todo.md`

## Where to start

1. Read `plan/00-overview.md` for full scientific rationale and study design
2. Check `plan/01-todo.md` for current tasks
3. Review `plan/02-decisions.md` for key analytical choices
4. Read `plan/03-analysis-protocol.md` for step-by-step analysis protocol
5. Read `plan/04-figure-plan.md` for planned manuscript figures

## Data sources

| Dataset | GEO | Cells | Description |
|---------|-----|-------|-------------|
| PC9 erlotinib time-series | GSE134839 | 6,508 (2,942 post-QC) | Primary: 6 time points (D0,D1,D2,D4,D9,D11) |
| PC9 osimertinib | GSE150949 | ~8,000 | Validation: 3rd-gen TKI |
| NSCLC mixed-lineage | GSE207422 | ~7,400 | Validation: patient tumors |

## Manuscript target

- **Primary:** Nature Communications / Genome Biology
- **Backup:** Cell Reports / Nucleic Acids Research
- **Preprint:** bioRxiv (submit simultaneously)
