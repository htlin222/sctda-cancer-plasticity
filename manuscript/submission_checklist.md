# Genome Biology Submission Checklist

**Target journal:** Genome Biology
**Article type:** Method
**Preprint:** bioRxiv (to deposit before journal submission)

---

## 1. Manuscript formatting requirements

### Genome Biology Method article structure
- [x] Title (concise, includes "framework" + biological scope)
- [x] Structured Abstract (Background / Results / Conclusions)
- [x] Keywords (6–10)
- [x] Background
- [x] Results (7 subsections)
- [x] Discussion
- [ ] Conclusions (ADD — currently merged with Discussion)
- [x] Methods
- [ ] Declarations block (CoI, funding, author contributions, ethics, data/code availability)
- [x] References (Vancouver style needed — currently numbered)
- [x] Figure captions
- [x] Table list

### Word count targets (Genome Biology Method)
- Background: 500–1,000 words (**currently ~500 — OK**)
- Results: 2,500–4,000 words (**currently ~2,500 — OK**)
- Discussion: 500–1,500 words (**currently ~800 — OK**)
- Methods: 1,500–3,000 words (**currently ~1,200 — EXPAND**)
- Total: 5,000–10,000 words

### Figures
- [x] Main figures: 5–8 (current: 8 figures — OK)
- [x] Resolution: 300 dpi vector PDF where possible
- [x] Font: Arial/Helvetica, sizes 7–9 pt
- [ ] Convert Unicode subscripts to LaTeX math for all captions (DONE in visualize.py, confirm exported figures)

### Supplementary
- [ ] Supplementary methods
- [ ] Supplementary figures (S1–S5 planned)
- [ ] Supplementary tables (S1–S4 planned)

---

## 2. Data availability statements (REQUIRED)

### GEO accessions (already public)
- GSE134839 (Aissa et al. 2021): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134839
- GSE150949 (Oren et al. 2021): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150949
- GSE243562 (Wurtz et al. 2024): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243562
- GSE131907 (Kim et al. 2020): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907

### Derived data (to deposit)
- [ ] Processed h5ad files → Zenodo (create record, get DOI)
  - pc9_erlotinib_analyzed.h5ad
  - pc9_osimertinib_processed.h5ad
  - pdx_ascl1_osimertinib_processed.h5ad
  - kim_luad_atlas_processed.h5ad
- [ ] Persistence diagrams (.npy files) → include in same Zenodo record
- [ ] Mapper graphs (.pkl files) → include in same Zenodo record

### Required data statement text
> All raw data are publicly available from the Gene Expression Omnibus under accessions GSE134839, GSE150949, GSE243562, and GSE131907. Processed AnnData files, persistence diagrams, and Mapper graphs generated in this study have been deposited at Zenodo (DOI: [pending]). The complete analysis pipeline, including all scripts and configuration files, is available at GitHub (https://github.com/<user>/sctda-cancer-plasticity) under the MIT license.

---

## 3. Code availability statement (REQUIRED)

### GitHub repository
- [ ] Create public GitHub repo
- [ ] Push all code (currently local at /Users/htlin/sctda-cancer-plasticity)
- [ ] Add README with installation and usage instructions
- [ ] Tag release v1.0.0 matching manuscript submission
- [ ] Register Zenodo integration to get DOI for v1.0.0

### Required code statement text
> All code is available at https://github.com/<user>/sctda-cancer-plasticity under the MIT license (version 1.0.0 archived at Zenodo DOI: [pending]). The repository includes all analysis scripts, the `sctda_plasticity` Python package, unit tests, and a Makefile orchestrating the full pipeline. The conda environment specification is in `envs/environment.yml`. Run `make pipeline` to reproduce all results from raw GEO data. Random seed 42 is used throughout.

---

## 4. Declarations block (REQUIRED by Genome Biology)

### Ethics approval
- [ ] N/A (all data are publicly available, no new human/animal samples collected)

### Consent for publication
- [ ] N/A

### Competing interests
- [ ] The authors declare that they have no competing interests.

### Funding
- [ ] [Add grant numbers, institutions]

### Author contributions (CRediT format)
- [ ] Conceptualization
- [ ] Data curation
- [ ] Formal analysis
- [ ] Investigation
- [ ] Methodology
- [ ] Software
- [ ] Supervision
- [ ] Validation
- [ ] Visualization
- [ ] Writing — original draft
- [ ] Writing — review & editing

### Acknowledgments
- [ ] Original data generators (Aissa et al., Oren et al., Wurtz et al., Kim et al.)
- [ ] Open-source software authors (scanpy, ripser, kmapper, persim, giotto-tda)
- [ ] Computing resources

---

## 5. Figures (current status)

| # | Title | File | Status |
|---|-------|------|--------|
| 1 | Study design + UMAP + EMT + pseudotime + PAGA | fig1_study_overview.pdf | ✅ Ready |
| 2a | Persistence barcodes by timepoint (GSE134839) | fig2_a_persistence_barcodes.pdf | ✅ Ready |
| 2b | H₁ persistence diagrams overlay | fig2_b_h1_persistence_diagrams.pdf | ✅ Ready |
| 2c | Permutation test histogram (D9) | fig2_c_permtest_D9.pdf | ✅ Ready |
| 2d | Wasserstein heatmap (discovery) | fig2_d_wasserstein_heatmap.pdf | ✅ Ready |
| 4a | Barcodes with CC genes | fig4_a_barcodes_with_cc.pdf | ⚠️ Need regenerate |
| 4b | Barcodes without CC genes | fig4_b_barcodes_without_cc.pdf | ⚠️ Need regenerate |
| 5 | Mapper + gene connectivity | (Mapper HTML + gene table) | ⚠️ Need composite |
| 6a | Cross-dataset H₁ diagram overlay | fig6_a_persistence_comparison.pdf | ✅ Ready |
| 6b | Osimertinib per-timepoint barcodes | fig6_b_validation_barcodes.pdf | ✅ Ready |
| 6c | Osimertinib internal Wasserstein | fig6_c_validation_wasserstein.pdf | ✅ Ready |
| 6d | H₁ over time both drugs | fig6_d_h1_over_time.pdf | ✅ Ready |
| 7a | PDX barcodes | fig7_a_pdx_barcodes.pdf | ✅ Ready |
| 7b | PDX H₁ diagrams | fig7_b_pdx_h1_diagrams.pdf | ✅ Ready |
| 7c | PDX Wasserstein heatmap | fig7_c_pdx_wasserstein.pdf | ✅ Ready |
| 7d | Kim patient atlas barcode | fig7_d_kim_atlas_barcodes.pdf | ✅ Ready |
| 8 | Master cross-system comparison | fig8_master_h1_comparison.pdf | ✅ Ready |

---

## 6. Tables (current status)

| # | Title | File | Status |
|---|-------|------|--------|
| 1 | Per-timepoint H₁ summary (GSE134839) | ph_summary_by_timepoint.csv | ✅ Ready |
| 2 | Validation permutation tests | permutation_tests_key_groups.csv | ✅ Ready |
| 3 | CC ablation comparison | cell_cycle_ablation_comparison.csv | ✅ Ready |
| S1 | CC ablation permtest | permutation_test_h1_ablation.csv | ✅ Ready |
| S2 | PDX gene connectivity | pdx_gene_connectivity_scores.csv | ✅ Ready |
| S3 | Mapper parameter sensitivity | mapper_sensitivity.csv | ✅ Ready |
| S4 | Full gene connectivity (cell line) | gene_connectivity_scores.csv | ✅ Ready |
| S5 | Wasserstein cross-dataset | wasserstein_validation_vs_discovery.csv | ✅ Ready |

---

## 7. To-do before submission

### Critical path (MUST do before submission)
1. [ ] Decide author list and affiliations
2. [ ] Create GitHub repo and push code
3. [ ] Get Zenodo DOI for code + data
4. [ ] Fill in Declarations block (CoI, funding, contributions)
5. [ ] Add a formal "Conclusions" section (split from Discussion)
6. [ ] Expand Methods to 1,500–3,000 words (add reproducibility subsection)
7. [ ] Regenerate Figure 4 (CC ablation barcodes with new font)
8. [ ] Compile composite Figure 5 (Mapper + gene attribution)
9. [ ] Convert references to Vancouver style
10. [ ] Final proofread by all authors

### Desirable but optional
- [ ] Add Figure 3 (dedicated Mapper graph visualization) — currently absorbed into Fig 5
- [ ] Supplementary methods document (computational details, parameter choices)
- [ ] Compile Word or LaTeX version of manuscript
- [ ] Internal review by 1–2 colleagues
- [ ] Preprint upload to bioRxiv 24–48 h before journal submission

### Preprint upload requirements (bioRxiv)
- [ ] Convert Markdown → PDF with pandoc
- [ ] Use standard bioRxiv LaTeX/Word template for camera-ready formatting
- [ ] Submit for deposit; acquire bioRxiv DOI
- [ ] Cite bioRxiv DOI in GB cover letter

---

## 8. Submission portal steps (Genome Biology)

1. Go to https://submission.springernature.com/
2. Select Genome Biology → Method article type
3. Upload:
   - Main manuscript (PDF + editable Word/LaTeX)
   - All figures (PDF/TIFF 300 dpi)
   - Tables (as separate file or embedded)
   - Supplementary materials (single PDF)
   - Cover letter (PDF)
4. Enter metadata:
   - Title, abstract, keywords
   - Author list, affiliations, ORCIDs, corresponding author
   - Conflict of interest, funding, ethical approval
   - Suggested and excluded reviewers
   - Data and code availability statements
5. Confirm preprint deposit status (yes, bioRxiv DOI)
6. Submit

---

## 9. Timeline estimate

| Task | Days |
|------|------|
| GitHub + Zenodo setup | 1 |
| Fill in declarations + authorship | 1 |
| Add Conclusions + expand Methods | 2 |
| Regenerate Figure 4 + composite Figure 5 | 1 |
| Convert references to Vancouver | 1 |
| Internal review | 3 |
| Final edits | 1 |
| Markdown → Word/LaTeX conversion | 1 |
| bioRxiv upload | 1 (wait for DOI) |
| Genome Biology submission | 1 |
| **Total** | **~2 weeks** |

---

## 10. Post-submission tracking

- [ ] Record submission date and manuscript ID
- [ ] Bookmark editorial tracking URL
- [ ] Initial decision expected: 2–6 weeks
- [ ] Major revision typical (Genome Biology): 2–4 months turnaround
- [ ] Plan for reviewer responses: set calendar reminder after 6 weeks
