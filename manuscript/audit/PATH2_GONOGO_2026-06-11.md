---
title: "Path 2 de-risking — GO/NO-GO decision"
date: 2026-06-11
decision: "NO-GO on Path 2 as specified (clonal state-entropy memory-decay)"
inputs: scripts/22-23, research memo (path2_research_memo.md), MeRLin GSE299589 inspection
tags: [go-no-go, path2, de-risk, decision]
---

# Path 2 de-risking — decision: NO-GO (as specified)

Proposed Path 2: a methods paper on a **clonal transcriptional-state memory statistic and its
decay under drug**, anchored on GSE150949 (PC9/osimertinib/Watermelon), validated externally,
with the cell-cycle-preserving lineage-shuffle null as the methodological hook.

Ran the full 4-step de-risk. **Three of four gates fail. Recommendation: NO-GO as specified;
fall back to Path 1 (the cautionary/methods paper).**

## Gate 1 — Novelty (research memo): MEDIUM-HIGH risk
Static "clones share state / within-clone heritability" is pre-empted (Lin et al. 2025/2026 Cell
Systems; PATH, Nat Genet 2024; Mathew 2023; MemorySeq). Only the *decay-under-drug trajectory +
cell-cycle-preserving null* is arguably novel — and narrowly. Survivable, but only if the other
gates hold. They don't.

## Gate 4 — Statistical robustness: FAILED (this is the killer)
The properly normalized statistic M(t)=1−H_obs/H_null is **not robust** on GSE150949:
- **Sign-flips:** decay D0→D14 > 0 in only **3/9** Leiden-resolution × clone-size settings. At
  k=2 (most data) decay is *negative* at every resolution (−0.505, −0.045, −0.166).
- **Wild instability:** M(D0) ranges −0.12 → 1.00 across settings; bootstrap 90% CIs cross zero
  and exclude the point estimates.
- **Non-monotonic trajectory:** M(D0)=0.48 → M(D3)=−0.02 → M(D7)=0.09 → M(D14)=0.25. A
  pre-registered monotonic-trend test (memo-required) fails.
- **Confounds in the trajectory:** D3 states track replicate (state~sample NMI=0.60 = batch
  artifact); D14 states partly track the high/med/low FACS sort (NMI=0.24).
- **Depth control halves it:** depth-stratified null cuts the D0→D14 decay 0.230 → 0.095.
- **Root cause (fundamental):** only 12–33 clones with ≥3 cells per timepoint, median clone
  size 3. You cannot estimate a within-clone state distribution from ~3 cells. GSE150949 is
  under-powered for this statistic, and no parameter choice fixes it.

This triggers the memo's explicit kill criterion ("depth/clone-size-matched null abolishes the
decay; confound, not memory").

## Gate 3 — External validation: BLOCKED (wrong modality)
The memo's #1 pick, **MeRLin / GSE299589, is bulk RNA-seq per clone** (STAR `ReadsPerGene.out.tab`
files, one bulk transcriptome per clone per timepoint; series type "Expression profiling by high
throughput sequencing", titled "[RNA-Seq]"). There is **no within-clone single-cell
distribution**, so the state-entropy statistic cannot be applied. The agent's 9/10 feasibility
rating did not verify modality. The MeRLin scRNA component is not in this series (only a
BioProject link, PRJNA1275025); locating a clone-rich scRNA sub-series is unverified archaeology.
Backups (GALILEO, SPLINTR) were not modality-verified either, so confidence is low.

## Gate 2 — Benchmarking: not reached (moot given gates 3,4).

## Why this is NO-GO, not "try harder"
The failure is not framing — it is **data power**. Lineage-traced scRNA with enough multi-cell
clones per timepoint AND a drug timecourse is rare, and the one strong candidate is bulk. The
in-house anchor itself sign-flips under standard robustness checks. Betting months on finding a
clone-rich external scRNA dataset that rescues a non-robust in-house signal is a bad wager.

## What IS still true and worth keeping
- Clonal state **memory exists at baseline** (M(D0)≈0.35–0.48, depth-matched, positive) — a real
  but unsurprising and pre-empted observation.
- The cell-cycle-preserving / dispersion-matched null framework is genuinely useful — but it is a
  **Path 1** contribution (a control toolkit), not a Path 2 discovery.

## Recommendation
1. **Path 1 (cautionary/methods paper)** is the robust, honest, publishable output. All its
   analyses are done (scripts 16–21). Fold in the Path 2 negative result as a second cautionary
   case: "naive clonal-memory statistics are also under-powered/confounded on standard
   lineage-tracing data; here are the controls."
2. **Revisit Path 2 only if** a genuinely clone-rich (many multi-cell clones/timepoint)
   lineage-traced scRNA drug-timecourse dataset is found and verified — a prerequisite, not a hope.
