---
title: "Consolidated Peer Review — TDA Cyclic Cell-State Plasticity Manuscript"
manuscript: "Topological data analysis reveals cyclic cell-state plasticity underlying drug tolerance in EGFR-mutant lung cancer"
date: 2026-06-10
review_type: simulated peer review (4-reviewer panel)
target_journal: Nature Communications
panel_verdict: Major revision
reviewers:
  - id: R1
    lens: Concept & Novelty
    recommendation: Major revision (leaning reject in current framing)
  - id: R2
    lens: Theory & Methodological Rigor
    recommendation: Major revision
  - id: R3
    lens: Data, Feasibility & Reproducibility
    recommendation: Major revision
  - id: R4
    lens: Editorial / Presentation / Synthesis
    recommendation: Major revision (consider transfer to Genome Biology / Cell Reports Methods)
files_reviewed:
  - manuscript/sections/abstract.tex
  - manuscript/sections/background.tex
  - manuscript/sections/results.tex
  - manuscript/sections/methods.tex
  - manuscript/sections/discussion.tex
  - manuscript/sections/conclusions.tex
  - manuscript/sections/figures_ncomms.tex
  - manuscript/sections/declarations.tex
  - manuscript/supplements/methods_extended.tex
  - manuscript/submission/checklist.md
  - manuscript/submission/declarations.md
tags: [peer-review, audit, manuscript, nature-communications, tda, persistent-homology]
---

# Consolidated Peer Review — "TDA reveals cyclic cell-state plasticity in drug-tolerant EGFR-mutant lung cancer"

**Panel verdict: Major revision (4/4 reviewers).** Two reviewers note it leans toward reject *in its current framing* but is salvageable. The reproducibility infrastructure and control battery are genuinely above field average; the problem is that the central claim is an interpretation laid on top of a fragile scalar, and several headline statements are contradicted by the paper's own numbers.

## The five issues every reviewer hit (these gate acceptance)

**1. Construct validity — a static loop is not temporal cycling.** This is the deepest problem (R1, R4). The biological claim is that cells *reversibly cycle* among states over time. The measurement is max H₁ persistence on a **single-timepoint snapshot point cloud** — a geometric property of the population's instantaneous distribution, not evidence that any cell traverses a loop, let alone reversibly. There is no circular-coordinate analysis, no persistent cohomology, no RNA-velocity-around-the-loop, and — most damningly — the osimertinib data carry **Watermelon lineage barcodes that could directly test cycling but are never used for it**. The same H₁ signal is fully consistent with a static branched/bimodal continuum. Either demonstrate true cycling (the lineage data is the obvious vehicle) or retreat the language from "cyclic plasticity" to "drug exposure increases the topological complexity of the state distribution."

**2. The headline statistic is fragile and unnormalized.** max H₁ = max(death − birth) is a single outlier-dominated bar, and the authors concede its scale is "metric-dependent" — yet they read raw magnitude increases as biology (R2, R3, R4). Their own numbers show factor-of-2–4 instability: CV up to 23.5% across just 5 seeds, factor-of-2 swing across PC counts, and the S/G2M regression "control" *increases* H₁ by up to **4.05×** (which they spin as confirmatory but reads as instability). No confidence intervals appear on any main-text value. Required: scale-normalized statistic as primary, bootstrap CIs on every reported number, and direction-significance at fixed cell count.

**3. The significance rests on a near-vacuous null; the meaningful one is underpowered.** The headline "p<0.01" comes from a gene-label null that shuffles genes *within each cell* and destroys all covariance — essentially any real data beats it, so it tests "is there any structure" not "is there a drug-dependent loop" (R2, R3). The one structure-preserving test (timepoint-label) was run at **n=50 permutations, for a single comparison**. This violates the project's own stated **n≥500** standard; validation cohorts used 100. Required: run the structure-preserving null at n≥500 for every headline comparison and replace the gene-label p-values.

**4. The heterogeneity confound + no negative control.** Treatment-naive *patient tumours* show the study's **highest** H₁ (4.77–5.37) — higher than fully drug-treated cell lines (3.78). If H₁ measured drug-induced cycling, untreated human tumours should not top the chart (R1). The parsimonious reading is that H₁ tracks population heterogeneity / number of coexisting states / cell count, all of which the authors concede grow with these variables. There is **no negative control** anywhere showing the statistic correctly reports *no* increase (R3). Required: dissociate H₁ from a heterogeneity/covariance-matched null and add a negative control.

**5. Overclaiming contradicted by the paper's own data.** The abstract says cyclicity "increases monotonically with drug exposure **across all systems**," but the erlotinib *discovery* cohort is **non-significant** (best p=0.098) and **non-monotonic** (D11=1.08 is the lowest of six timepoints) (R1, R4). The abstract's "replicates in a 14-patient longitudinal cohort" oversells a pooled trend plus **n=3** matched biopsies that the Discussion itself calls "hypothesis-generating" (R3, R4). Fix the abstract/title to match the evidence (monotonic under osimertinib in three systems; underpowered under erlotinib).

## Additional substantive findings

- **Internal numerical inconsistencies (data-integrity flag).** The YU-006 PDX untreated baseline appears as **four different values** — 2.36 / 2.81 / 2.89 / 3.62 — across results/methods/tables. Osimertinib D0 is **1.41 in Results but 1.71 in the Harmony section**. A reader cannot trace any number to a definitive source. Needs one master table with cell count + seed per value (R2, R3).
- **EMT circularity.** The cell-line Mapper uses the **EMT score as its filter**, then reports EMT genes as top-connected and Hallmark-EMT as top-enriched — partly tautological. The unbiased PC1 filter shows cell-cycle markers co-ranking, softening the clean "EMT anchors the loop" story (R1, R2).
- **PDX mouse-gene removal** strips 29% of features by a lowercase-symbol regex the authors admit is inferior to alignment-based methods — bidirectional misassignment directly perturbs the PCA the PDX results (the key *in vivo* validation) depend on (R1, R2, R3).
- **"Bit-stable across runs"** is contradicted by the paper's own S14 (BLAS/OpenMP nondeterminism) and the 23.5% CV. **Version tag** says v2.7.1-submission but repo tags reach v2.7.7 (R3).
- **The titular "loop" is never shown as a cycle** — no Mapper graph containing a demonstrated loop is displayed; the planned Mapper-graph figure (Fig 3) and circular-profile figure were **dropped**, and those were exactly the panels that would substantiate the claim (R4).

## Presentation / editorial

- **Not submittable to Nature Communications as-is:** 7,835 words vs 5,000 limit; 14 display items vs 10; community `nature.cls` instead of official `sn-jnl.cls`; funding placeholder unresolved; Zenodo DOIs "to follow" (NC requires an archived DOI, not a GitHub tag).
- **AI-image declaration contradiction:** `declarations.md`, `checklist.md`, and the Fig 1 caption disagree on whether Figure 1a contains a gpt-image-1 render — reconcile before submission (integrity flag).
- **Missing literature:** Shaffer 2017 (foundational non-genetic resistance) is in the .bib but uncited; CellRank/scVelo (the obvious "why not just use this" for cyclic dynamics) not distinguished; persistence-image/landscape vectorized statistics (Adams 2017, Bubenik 2015) — directly relevant to the single-scalar weakness — not cited.

## Venue fit (panel consensus)

The defensible contribution is a **reproducible, well-controlled framework/resource** for TDA on perturbation scRNA-seq — not a new mechanism (the cycling biology was already established by the lineage-tracing papers they cite). As a biological-discovery paper for Nature Communications the novelty is thin and the statistic fragile; as a methods/resource paper it is strong. Two reviewers recommend transferring to **Genome Biology** (a cover letter for it already exists in the repo) or **Cell Reports Methods**. NC is defensible only if issues 1 and 2 are decisively answered.

## Highest-impact next move

If keeping the NC framing: use the **Watermelon lineage barcodes in GSE150949** to show the same lineages occupy different angular positions over time. That single analysis would convert "static correlate" into "demonstrated cycling" and simultaneously answer the construct-validity, loop-is-real, and heterogeneity objections.

## Per-reviewer recommendations

| Reviewer | Lens | Recommendation |
|----------|------|----------------|
| R1 | Concept & Novelty | Major revision (leaning reject in current framing) |
| R2 | Theory & Methodological Rigor | Major revision |
| R3 | Data, Feasibility & Reproducibility | Major revision |
| R4 | Editorial / Presentation / Synthesis | Major revision (consider transfer to Genome Biology / Cell Reports Methods) |
