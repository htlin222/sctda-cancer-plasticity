# Cover Letter — Genome Biology submission

**To:** The Editors, *Genome Biology*
**From:** Hsieh-Ting Lin, M.D. — Department of Oncology, Koo Foundation Sun Yat-Sen Cancer Center, Taipei, Taiwan (htlin222@kfsyscc.org)
**Subject:** Manuscript submission — "Drug tolerance takes shape: persistent homology reveals cyclic cell-state plasticity in EGFR-mutant lung cancer"

---

Dear Editors,

I am pleased to submit for your consideration the manuscript **"Drug tolerance takes shape: persistent homology reveals cyclic cell-state plasticity in EGFR-mutant lung cancer"** as a Method article in *Genome Biology*.

## Why this work fits *Genome Biology*

Standard single-cell RNA-seq trajectory inference methods — UMAP, Monocle, Slingshot, PAGA, RNA velocity — assume that cell states lie along tree-like or directed manifolds. This assumption fails for reversibly plastic states such as those observed in drug-tolerant persister (DTP) cells, EMT, and other stress-response phenotypes, where cyclic or multi-connected topology is expected. Existing TDA work on scRNA-seq (Rizvi et al. 2017) demonstrated feasibility in a developmental context but did not provide an open, validated pipeline with formal null-model testing or cross-dataset benchmarking.

Our work fills this gap. We present an **open-source, modular framework** that combines persistent homology (H₀, H₁), the Mapper algorithm, permutation-based null tests, cell-cycle ablation controls, Wasserstein distance comparisons, and Mapper-based gene connectivity attribution into a single reproducible pipeline. The framework is applied without modification across four heterogeneous scRNA-seq datasets spanning **cell lines, patient-derived xenografts, and patient tumor atlases**, covering 258,355 cells pre-QC and demonstrating that a single topological statistic (max H₁ persistence) reproducibly captures drug-induced cell-state reorganization. All code is lint-clean, test-covered (8/8 tests pass), and published under MIT license.

We believe the manuscript matches *Genome Biology*'s mission on three specific criteria:

1. **Methodological novelty with rigorous reproducibility.** The pipeline is fully containerized via a single `make pipeline` target, uses deterministic seeding, includes unit tests, and has been validated on four independent publicly available GEO datasets. We specifically emphasize negative controls (cell-cycle ablation, permutation null) and honest limitations discussion.

2. **Cross-scale validation.** Our results span cell line, PDX, and patient tumor baselines, demonstrating the framework's robustness across biological systems without assuming equivalence between them. This is a rare level of validation for a methods paper.

3. **Immediate reusability.** The framework accepts standard `AnnData` input and exposes each analytical stage as an independent module. Users can substitute their own datasets, swap filter functions, or adjust parameters without touching pipeline internals. We anticipate uptake in drug-screening, developmental, and immunology contexts.

## Significance and claims

### Central methodological claim
A scale-independent topological statistic (max H₁ persistence on Vietoris-Rips filtration of top-30 PCA space) detects and quantifies drug-induced cell-state plasticity reproducibly across independent scRNA-seq datasets.

### Central biological claim
In four independent EGFR-mutant NSCLC datasets, maximum H₁ persistence increases monotonically with osimertinib exposure (2.7× over 14 days in PC9 cell line; 1.4× to 2.2× in residual-disease PDX models), consistent with the reorganization of cell-state space toward the baseline topological complexity observed in treatment-naive patient tumors. Loop-associated genes are enriched for EMT programs in a convergent but system-specific manner.

### Independent evidence for cyclic DTP plasticity
Our topological signature independently validates the cycling-persister phenomenon recently reported by Oren et al. (2021, *Nature*) using Watermelon lineage tracing, while extending the observation to two independent PDX systems and quantifying it with a reproducible statistic.

## Why we chose *Genome Biology* over other venues

*Genome Biology* has a strong history of publishing methods papers that combine mathematical rigor with broad biological applicability (e.g., Rizvi et al. 2017 scTDA, Saelens et al. 2019 dynverse trajectory benchmark, Lopez et al. 2018 scVI). Our framework is designed to be a similar contribution — a pipeline the community can pick up, apply, and extend, rather than a single-cohort finding paper. The readership of *Genome Biology*, which includes both method developers and biologists who apply methods, is the ideal audience for a cross-validated TDA framework.

## Suggested reviewers (optional)

We suggest reviewers with expertise in TDA, scRNA-seq trajectory inference, and lung cancer drug tolerance. Please let us know if you would like us to provide specific suggestions.

## Conflict of interest and dual submission

The authors declare no conflict of interest. This manuscript has been deposited on bioRxiv as a preprint [DOI pending upon submission] and is not under consideration at any other journal.

We look forward to your assessment and are happy to provide any additional information the editorial team or reviewers may request.

Yours sincerely,
[Corresponding author signature block]
