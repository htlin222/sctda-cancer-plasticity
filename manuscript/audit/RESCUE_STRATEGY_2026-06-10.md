---
title: "Publication rescue strategy — from 'cyclic plasticity' to a dispersion-controlled topological signature"
date: 2026-06-10
status: empirical core locked; manuscript reframe pending
analyses: scripts/16_lineage_circular_coords.py, scripts/17_dispersion_null_and_stability.py
evidence: results/foundation/
tags: [rescue, reframe, strategy, manuscript, tda]
---

# Publication rescue strategy

## What died and what survived

The original thesis — drug tolerance is **reversible cyclic plasticity** detected as rising
max-H₁ — has two fatal problems confirmed by direct analysis:

1. **Cells do not traverse the loop** (construct validity). Validated circular coordinates
   (dreimac DSPVJ) collapse to one angle in all 5 cohorts (Rayleigh R ≥ 0.965 vs 0.881 for a
   12%-occupied control loop); Watermelon clones are *tighter* than random, not spread.
   → "cycling" is unsupported.
2. **At baseline, max-H₁ is just dispersion.** Against a covariance-matched Gaussian null
   (same mean/cov, no real structure), baseline ratios are ≈0.9 (real ≤ null).

**What survives — and is genuinely publishable:**

> Under EGFR-TKI (osimertinib) in cell-line and PDX models, drug-tolerant persister emergence
> is accompanied by **genuine non-Gaussian topological reorganization** of the cell-state space:
> max-H₁ rises **above a dispersion-matched null specifically under drug** (not at baseline),
> robustly across subsamples. The signal is a **sparse persistent loop occupied by a minority
> of cells**, not population-wide cycling.

## Locked evidence (scripts/17, bootstrap 30× / Gaussian null)

| group | real H₁ [90% CI] | gauss | ratio | real>gauss | reading |
|---|---|---:|---:|---:|---|
| osi D0 | 1.41 [1.31,1.60] | 1.49 | 0.95 | 31% | baseline = dispersion |
| osi D3 | 1.82 [1.82,2.15] | 2.00 | 0.91 | 18% | dispersion |
| osi D7 | 2.67 [2.45,3.04] | 2.04 | 1.31 | **100%** | structure emerges |
| osi D14 | 3.78 [2.99,3.78] | 2.66 | 1.42 | **98%** | structure |
| erl D0/D9/D11 | 1.49 / 1.77 / 1.08 | ~1.6 | 0.93/0.99/0.70 | 33/42/0% | **no excess — drop** |
| PDX untreated | 2.79 [2.23,4.03] | 2.30 | 1.21 | 86% | marginal |
| PDX residual | 5.46 [3.58,7.34] | 2.80 | 1.95 | **99%** | strong drug effect |
| Maynard TN/PER/PD | 5.49 / 6.04 / 5.58 | ~3.6 | 1.52/1.69/1.51 | 99/98/95% | excess at ALL stages incl naive → heterogeneity, not drug |

Directional robustness: osi D7>D0 and D14>D0 in **100%** of bootstrap pairs; PDX residual>untreated **97%**.

## Keep / drop / reframe

- **KEEP:** osimertinib cell-line + PDX drug contrasts; the Gaussian-matched null as the
  central method; bootstrap directional robustness; cell-cycle controls (still valid);
  reproducible pipeline.
- **DROP / DEMOTE:** "cyclic / reversible plasticity" language and the donut-cycling concept
  figure; the erlotinib "discovery" as a positive result (it shows no excess structure — recast
  as an underpowered/negative cohort, which actually strengthens specificity); cross-system
  "moves toward patient topology" master claim.
- **REFRAME:** patient (Maynard/Kim) data as showing that *baseline tumor heterogeneity*
  already produces excess topological structure → the clean drug-induction signal is the
  controlled in-vitro/PDX systems. Honest, and answers reviewer #1's "naive tumors are highest".

## How this answers each reviewer objection

- **#1 construct validity (cycling):** removed; claim is now "topological reorganization,"
  explicitly distinguished from cycling, with the circular-coordinate analysis shown as a
  negative control proving we are NOT overclaiming traversal.
- **#2 dispersion confound + vacuous null:** the covariance-matched Gaussian null IS the fix and
  the headline method; gene-label null retired.
- **#2 outlier-driven / unnormalized:** trim-1% test (7/9 robust) + bootstrap CIs reported.
- **#3 instability:** point estimates CV 6–45% disclosed; claims rest on directional bootstrap
  robustness, not point values.
- **#3 numerical inconsistencies:** rebuild one master table from scripts/17 output (single seed/source).
- **#1/#4 overclaim:** scope narrowed to osi+PDX; erlotinib/patient honestly bounded.

## Remaining work

1. **Biology** (in progress): recover gene symbols (processed h5ad lost them → use raw CSV),
   identify the minority cells driving the non-Gaussian structure in osi D14 / PDX residual,
   test for a coherent drug-tolerant/EMT/quiescence program controlling for count depth.
2. **Manuscript reframe** (orchestrated): abstract, title, intro, results, discussion, figures,
   reviewer response — all grounded in the locked table above.
3. **Negative control** is now built-in (erlotinib + baseline + Gaussian null all serve as nulls).
4. **Provenance fix:** processed objects lost gene names — note and regenerate with symbols.
