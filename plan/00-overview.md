# Project Plan: scTDA Cancer Plasticity

Generated: 2026-04-14
Status: Approved
Methodology: Data Science / Computational Biology — Publication Pipeline

---

## Central Hypothesis

Standard trajectory inference methods (Monocle, Slingshot, PAGA, RNA velocity) assume
tree-like or directed structures in cell-state space. However, drug-tolerant cancer cells
exhibit phenotypic plasticity — reversible transitions between states (e.g., EMT ↔ MET,
epithelial ↔ mesenchymal) — that would manifest as **loops (H₁ features)** in the
expression topology. These loops are invisible to tree-based methods but detectable by
persistent homology.

**If we can show:**
1. Statistically significant H₁ loops exist in drug-tolerant cell populations but not in untreated cells
2. These loops are NOT merely cell-cycle artifacts
3. The loops correspond to biologically meaningful reversible state transitions (EMT/MET, metabolic switching)
4. The finding replicates across independent datasets and drug types

**Then we have a high-impact paper** demonstrating that drug tolerance operates through cyclic plasticity, not linear fate transitions.

---

## Study Design

### Discovery Cohort: PC9 + Erlotinib (GSE134839)

EGFR-mutant NSCLC cell line PC9 treated with erlotinib (2 μM) at multiple time points.
Drop-seq data from Aissa et al., Nature Communications 2021.

- Day 0 (untreated): baseline topology
- Day 1: early response
- Day 3: drug-tolerant persisters emerging
- Day 11: expanded persisters (DTEPs)
- Drug holiday: cells removed from drug after treatment

This longitudinal design is ideal for tracking topological evolution.

### Validation Cohort 1: PC9 + Osimertinib (GSE150949)

Same cell line, different TKI (3rd-generation). If the same H₁ loop appears,
it's a general drug tolerance phenomenon, not erlotinib-specific.

### Validation Cohort 2: Patient NSCLC Mixed-Lineage (GSE207422)

Patient-derived NSCLC with mixed-lineage features (ADC + SCC markers in same cells).
If patient tumor cells show similar topological signatures, the finding has clinical relevance.

---

## Analysis Architecture

```
[Raw Data] ──→ [QC + Preprocessing] ──→ [Standard scRNA-seq]
                                              │
                                   ┌──────────┴──────────┐
                                   │                      │
                           [Baseline Analysis]    [Topological Analysis]
                           (UMAP, clustering,     (PH, Mapper, gene
                            trajectory, velocity)  connectivity)
                                   │                      │
                                   └──────────┬──────────┘
                                              │
                                   [Comparative Analysis]
                                   (What does TDA see that
                                    standard methods miss?)
                                              │
                                   [Biological Validation]
                                   (Gene attribution, pathway
                                    enrichment, cell cycle ablation)
                                              │
                                   [Cross-Dataset Validation]
                                              │
                                   [Manuscript Figures & Text]
```

---

## Phase 1 Features (Minimum Publishable Unit)

- [x] Project scaffold and plan
- [ ] Data download pipeline (GSE134839, GSE150949)
- [ ] Preprocessing pipeline (QC → normalize → HVG → PCA)
- [ ] Standard baseline analysis (UMAP, Leiden clustering, marker genes)
- [ ] Persistent homology on PCA space (H₀, H₁, H₂) per time point
- [ ] Permutation test for H₁ significance
- [ ] Mapper visualization with multiple filter functions
- [ ] Cell cycle ablation test
- [ ] Gene connectivity analysis on significant topological features
- [ ] Pathway enrichment (GO, Hallmark gene sets)
- [ ] Comparison figure: what TDA sees vs. what UMAP/velocity miss
- [ ] Validation on GSE150949

## Phase 2 (Strengthening)

- [ ] Persistent cohomology for circular coordinates on loops
- [ ] Patient tumor validation (GSE207422)
- [ ] Multi-distance-metric robustness analysis
- [ ] Subsampling stability analysis
- [ ] Integration with drug sensitivity assay data if available

## Phase 3 (Manuscript)

- [ ] Main figures (6-7 panels)
- [ ] Supplementary figures (15-20)
- [ ] Methods section with full reproducibility details
- [ ] Draft manuscript in LaTeX or Markdown
- [ ] bioRxiv preprint

---

## Chosen Stack

| Component | Choice | Rationale |
|-----------|--------|-----------|
| Language | Python 3.11 | scanpy ecosystem |
| Environment | conda/mamba | scientific reproducibility |
| scRNA-seq | scanpy 1.10+ | standard, well-documented |
| Persistent homology | ripser 0.6+ | fastest PH implementation |
| PH visualization | persim | persistence diagrams, landscapes |
| Mapper | kepler-mapper 2.0+ | scikit-learn compatible |
| TDA framework | giotto-tda 0.6+ | end-to-end TDA+ML |
| Enrichment | gseapy, decoupler | pathway analysis |
| Velocity | scvelo | RNA velocity comparison |
| Figures | matplotlib + seaborn | publication quality |

---

## Risk Register

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|------------|
| No significant H₁ in any dataset | Medium | Fatal | Screen 3+ datasets before committing |
| H₁ is only cell cycle | High | Major | Cell cycle ablation is a core analysis step |
| Reviewer says "just noise" | Medium | Major | Permutation tests, multiple filtrations, cross-validation |
| Computational bottleneck | Low | Minor | Subsample for exploration, full data for final results |
| Scooped | Low | Major | Preprint on bioRxiv as soon as core result is solid |

---

## Timeline

| Week | Phase | Milestone |
|------|-------|-----------|
| 1-2 | Setup | Environment, data download, preprocessing |
| 3-4 | Baseline | Standard scRNA-seq analysis complete |
| 5-7 | Core TDA | PH + Mapper + statistical tests |
| 8-9 | Biology | Gene attribution + pathway enrichment |
| 10-11 | Validation | Second dataset confirms findings |
| 12-14 | Figures | Publication-quality figures |
| 15-17 | Writing | Manuscript draft |
| 18 | Submit | bioRxiv + journal submission |

---

## Open Questions (to resolve during analysis)

| Question | Status | Decision |
|----------|--------|----------|
| Which distance metric for PH? | Open | Test Euclidean, correlation, diffusion; pick most robust |
| How many PCs to use? | Open | Elbow plot + test 20/30/50; check PH stability |
| Subsample size for PH? | Open | Start with all cells if <5000; subsample + bootstrap if larger |
| Which Mapper filter function? | Open | Test PC1, diffusion component 1, EMT score |
| Cell cycle gene list source? | Open | Tirosh et al. 2016 (standard) |
