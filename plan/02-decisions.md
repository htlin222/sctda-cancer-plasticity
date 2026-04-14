# Architecture Decision Log

## ADR-001: Use PCA space for persistent homology, not raw expression

- **Date:** 2026-04-14
- **Status:** Accepted
- **Context:** Persistent homology on raw expression (~20,000 genes) is computationally
  infeasible and noisy. Need a reduced representation.
- **Decision:** Compute PH on the top 30-50 PCA components of the HVG-filtered,
  log-normalized expression matrix. Test stability across 20, 30, and 50 PCs.
- **Consequences:** PCA is linear, so we may miss nonlinear structure. However, this is
  standard in the field and makes results comparable. We mitigate by also testing PH
  on diffusion map coordinates (nonlinear) as a robustness check.

## ADR-002: Permutation strategy for H₁ significance

- **Date:** 2026-04-14
- **Status:** Accepted
- **Context:** Need to distinguish real H₁ loops from noise. Following scTDA (Rizvi et al. 2017)
  approach of permuting gene labels.
- **Decision:** For each of 500 permutations, independently permute gene labels for each cell,
  rebuild PCA, compute H₁. The p-value for an observed H₁ feature is the fraction of
  permutations producing a feature with death time ≥ observed.
- **Consequences:** Conservative test. 500 permutations gives resolution to p=0.002.
  Computationally feasible on Mac Mini for <5000 cells.

## ADR-003: Cell cycle gene list

- **Date:** 2026-04-14
- **Status:** Accepted
- **Context:** Cell cycle is the most common confound producing H₁ loops. Need a definitive
  gene list for ablation experiments.
- **Decision:** Use Tirosh et al. 2016 (Science) S-phase and G2/M-phase gene lists, which
  are the standard in scanpy (sc.tl.score_genes_cell_cycle). Total ~97 genes.
- **Consequences:** Well-established, widely accepted by reviewers. May miss some
  cell-cycle-related genes, but conservative is better for our argument.

## ADR-004: Primary dataset choice

- **Date:** 2026-04-14
- **Status:** Accepted
- **Context:** Need a dataset with (a) time series, (b) known plasticity, (c) manageable size,
  (d) published baseline for comparison.
- **Decision:** GSE134839 (PC9 + erlotinib, Aissa et al. Nat Commun 2021). Time-series
  Drop-seq with drug-tolerant persister emergence. Already has RNA velocity analysis
  as baseline comparison. ~5000 cells total.
- **Consequences:** Cell line data limits clinical generalizability. We address this with
  patient tumor validation dataset in Phase 2. Cell line advantage: low noise, clear
  experimental design.

## ADR-005: Mapper filter function strategy

- **Date:** 2026-04-14
- **Status:** Accepted
- **Context:** Mapper output depends heavily on filter function choice. Need principled approach.
- **Decision:** Run Mapper with three filter functions:
  1. PC1 (captures largest variance axis — standard choice)
  2. EMT score (biologically motivated — tests EMT loop hypothesis directly)
  3. Diffusion component 1 (nonlinear, captures transition dynamics)
  Report all three; focus narrative on the one that best reveals the loop.
- **Consequences:** Multiple comparisons concern, but this is exploratory visualization,
  not hypothesis testing. PH provides the statistical backbone.

## ADR-006: Distance metric for persistent homology

- **Date:** 2026-04-14
- **Status:** Proposed
- **Context:** Euclidean distance in PCA space is standard but may not capture gene
  expression similarity well. Correlation distance is scale-free.
- **Decision:** Primary analysis uses Euclidean in PCA space (most comparable to literature).
  Supplementary analysis uses correlation distance for robustness.
- **Consequences:** If results agree across metrics, strengthens the claim.
