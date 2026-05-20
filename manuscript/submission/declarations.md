# Declarations

## Competing interests

The authors declare no competing financial or non-financial interests.

## Ethics

This study analyses **publicly available, de-identified** single-cell RNA-seq data and does not involve new human-subject or animal experiments.

| Dataset | GEO / source | Original IRB / approval |
|---------|--------------|-------------------------|
| PC9 erlotinib time-series | GSE134839 (Aissa et al., *Nat Commun* 2021) | Cell-line work; no IRB needed |
| PC9 osimertinib (Watermelon) | GSE150949 (Oren et al., *Nature* 2021) | Cell-line work; no IRB needed |
| EGFR PDX osimertinib | GSE243562 (Wurtz et al. 2024) | Original study IRB-approved (mouse PDX) |
| LUAD atlas (Kim 2020) | GSE131907 | Original study IRB-approved (Kim et al., *Nat Commun* 2020) |
| EGFR-mutant longitudinal cohort | Maynard et al. *Cell* 2020 (public archive) | Original study IRB-approved at UCSF |

Re-analysis of these public datasets is consistent with the data-use terms of each contributing study.

## Data availability

- **Raw scRNA-seq data:** GEO accessions GSE134839, GSE150949, GSE243562, GSE131907; Maynard *et al.* 2020 EGFR-mutant cohort via the authors' public Google Drive archive cited in that paper.
- **Processed AnnData objects, persistence diagrams, and Mapper graphs:** regenerable end-to-end via `make pipeline` (see Code availability). A version-tagged release of intermediate outputs will be deposited at **Zenodo** prior to acceptance and the DOI added here.
  - *Zenodo DOI:* `<to be assigned at submission>`.
- **Supplementary Tables S1–S11** are provided as separate files at submission.

## Code availability

- **Repository:** <https://github.com/htlin222/sctda-cancer-plasticity> (MIT license).
- **Reproducibility:** the complete pipeline runs from raw GEO data to all figures and tables via `make pipeline`. A `CITATION.cff` provides machine-readable software-citation metadata.
- **Environment:** conda environment in `envs/environment.yml`; tested on macOS arm64 and Linux x86_64.
- **Versioned release:** a release tag following the existing `vX.Y.Z-submission` convention (current: `v2.7.0-submission` plus the NC-prep delta — final tag TBD) and a Zenodo software DOI will be created at submission (`<DOI to be assigned>`).

## Funding

This work received no specific external grant funding and was carried out using the authors' institutional computing resources at Koo Foundation Sun Yat-Sen Cancer Center (H.-T.L.) and the Center for Computational and Systems Biology, National Taiwan University (Y.-H.T.).

## AI-tool usage statement

- Software-development assistance: **Claude Code** (Anthropic) was used for code-pipeline scaffolding, refactoring, and figure-script implementation. All scientific decisions, statistical analyses, and prose are the authors' own.
- **Figure 1 panel (a)** is a photorealistic 3D-render-style composite illustrating topological equivalence (coffee mug ≅ donut-with-handle ≅ donut). The image was generated with **OpenAI `gpt-image-1`** on a chroma-key background, matted to alpha, and embedded into the figure via matplotlib (`scripts/fig1_concept.py`). The image carries no data and is illustrative only. We acknowledge that *Nature Communications* generally restricts AI-generated imagery; we use this asset to communicate a purely conceptual topological-equivalence relationship for the broad non-specialist audience and will replace it with a hand-illustrated or vector schematic on editorial request. Panels (b) and (c) are programmatically generated via matplotlib from non-AI data (`scripts/fig1_concept.py`).

## Inclusion & ethics statement

This research did not involve human participants or vulnerable populations. The five publicly available datasets re-analysed here originate from studies conducted in compliance with each study's institutional ethics approval. The two authors are based in Taiwan; no exclusionary collaboration practices apply.
