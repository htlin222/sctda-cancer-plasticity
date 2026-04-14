# scTDA Cancer Plasticity

> **Drug tolerance takes shape: persistent homology reveals cyclic cell-state plasticity in EGFR-mutant lung cancer.**
> An open-source topological framework validated across cell lines, patient-derived xenografts, and patient tumors.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/)
[![Preprint PDF](https://img.shields.io/badge/preprint-PDF-red)](manuscript/preprint.pdf)

**Preprint:** [`manuscript/preprint.pdf`](manuscript/preprint.pdf) (bioRxiv DOI: pending)
**Target journal:** *Genome Biology* (Method article)

---

## What this is

An open-source pipeline that combines persistent homology ($H_0$, $H_1$) and the Mapper algorithm with cell-cycle ablation, permutation null testing, Wasserstein diagram comparison, and Mapper-based gene attribution into a single reproducible framework for longitudinal scRNA-seq analysis.

Applied to four public EGFR-mutant NSCLC datasets covering erlotinib, osimertinib, PDX models, and a patient LUAD atlas.

## Headline result

The maximum $H_1$ persistence — a scale-independent topological statistic — increases monotonically with osimertinib exposure across cell line, PDX, and patient scales, converging toward the baseline topological complexity of treatment-naive patient tumors. Loop-associated genes are enriched for EMT programs. See [`manuscript/preprint.pdf`](manuscript/preprint.pdf).

## Quick start

```bash
git clone git@github.com:<user>/sctda-cancer-plasticity.git
cd sctda-cancer-plasticity

# 1. Install environment (Python 3.11, ~15 min)
mamba env create -f envs/environment.yml
conda activate sctda
pip install -e .

# 2. Download data (public GEO; ~1 GB)
make download

# 3. Run the pipeline (cell line + osimertinib validation; ~30 min)
make pipeline

# 4. Extended validation (PDX + Kim atlas; ~20 min)
python scripts/08_validation_pdx.py
python scripts/09_validation_kim_atlas.py
python scripts/10_ablation_permtest.py
python scripts/11_pdx_gene_attribution.py

# 5. Build the manuscript PDF
cd manuscript && make pdf
```

## Pipeline stages

| # | Script | Purpose |
|---|--------|---------|
| 01 | `scripts/01_preprocess.py` | QC, normalize, HVG, PCA |
| 02 | `scripts/02_standard_analysis.py` | UMAP, Leiden, EMT score, DPT, PAGA |
| 03 | `scripts/03_topology.py` | Per-group PH, Mapper, Wasserstein, permutation tests |
| 04 | `scripts/04_cell_cycle_ablation.py` | $H_1$ with/without cell-cycle genes (ratio test) |
| 05 | `scripts/05_gene_attribution.py` | Mapper gene connectivity, DE, enrichment |
| 06 | `scripts/06_validation.py` | GSE150949 osimertinib cohort |
| 07 | `scripts/07_figures.py` | Composite publication figures |
| 08 | `scripts/08_validation_pdx.py` | ASCL1 PDX (GSE243562) $H_1$ comparison |
| 09 | `scripts/09_validation_kim_atlas.py` | Kim 2020 patient LUAD baseline (GSE131907) |
| 10 | `scripts/10_ablation_permtest.py` | Permutation null on CC-ablated data |
| 11 | `scripts/11_pdx_gene_attribution.py` | Cross-system gene-overlap analysis |
| 12 | `scripts/12_convert_tables_to_parquet.py` | CSV → parquet for large tables |

## Datasets

| Dataset | GEO | Cells (post-QC) | Role |
|---------|-----|-----------------|------|
| Aissa 2021 PC9 + erlotinib | [GSE134839](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134839) | 2,942 | Discovery (6 timepoints) |
| Oren 2021 PC9 + osimertinib | [GSE150949](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150949) | 4,999 (subsampled) | Validation (4 timepoints) |
| Hu 2024 EGFR-mutant PDX | [GSE243562](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243562) | 18,874 | PDX validation (pre/post) |
| Kim 2020 patient LUAD atlas | [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907) | 4,513 (subsampled) | Treatment-naive reference |

## Repository structure

```
.
├── manuscript/              # LaTeX manuscript (source of truth)
│   ├── preprint.tex         # Main document
│   ├── references.bib       # Validated citations
│   ├── preprint.pdf         # Compiled (22 pages)
│   ├── cover_letter_genome_biology.md
│   ├── submission_checklist.md
│   └── Makefile             # `make pdf` / `make biorxiv`
├── src/sctda_plasticity/    # Core Python package
│   ├── data.py              # Dataset loaders and QC
│   ├── topology.py          # PH, Mapper, gene connectivity
│   ├── visualize.py         # Publication-quality plotting
│   └── config.py            # Parameters and paths
├── scripts/                 # Pipeline stages (01-12)
├── tests/                   # 8 pytest unit tests (TDA synthetics)
├── envs/environment.yml     # Conda env (Python 3.11 + scanpy + ripser)
├── results/
│   ├── figures/             # 15+ publication figures (PDF + PNG)
│   ├── tables/              # Summary tables (parquet + CSV)
│   ├── ph/                  # Persistence diagrams (.npy)
│   ├── mapper/              # Mapper graphs (.pkl + HTML)
│   └── enrichment/          # gseapy pathway outputs
├── plan/                    # Project plan & decisions
├── Makefile                 # Top-level orchestration
├── LICENSE                  # MIT
└── README.md                # This file
```

## Data availability

| Tier | Where | How to get it |
|------|-------|---------------|
| Raw GEO files | public GEO | `make download` (~1 GB) |
| Processed AnnData h5ad | Zenodo (DOI pending) | `wget <url>` or regenerate with `make pipeline` |
| Persistence diagrams, Mapper graphs, figures, small tables | This repo | already included |
| Large tables (gene connectivity, DE) | This repo as parquet | `pd.read_parquet(...)` |

Parquet files are 42-66% smaller than CSV; original CSVs are regenerable. See `.gitignore` for details.

## Reproducibility

- **Environment**: Python 3.11 pinned in `envs/environment.yml`
- **Seeds**: all stochastic ops use `random_state=42`
- **Tests**: `pytest` (8/8 pass on macOS arm64 and Linux x86_64)
- **Lint**: `ruff check` clean (100-char line limit)
- **One command**: `make pipeline` reproduces all figures and tables from raw GEO data

## Citation

If you use this framework, please cite the preprint:

```bibtex
@article{sctda_plasticity_2026,
  title   = {Drug tolerance takes shape: persistent homology reveals cyclic cell-state plasticity in EGFR-mutant lung cancer},
  author  = {...},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {pending}
}
```

Also cite the original data generators: Aissa et al. 2021 (PMID:33712615), Oren et al. 2021 (PMID:34381210), Hu et al. 2024 (PMID:38359163), Kim et al. 2020 (PMID:32385277).

## License

MIT. See [LICENSE](LICENSE).

## Contact

Pull requests and issues welcome. Corresponding author: [email pending]
