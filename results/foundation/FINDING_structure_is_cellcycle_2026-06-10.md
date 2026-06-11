---
title: "Decisive finding — the drug-induced topological 'structure beyond dispersion' is largely the cell-cycle loop"
date: 2026-06-10
analyses: scripts/18 (biology), scripts/19 (per-tp control), scripts/21 (pooled CC-regression, decisive)
dataset: GSE150949 PC9 osimertinib (primary system)
status: supersedes the optimistic rescue framing in RESCUE_STRATEGY_2026-06-10.md
tags: [foundation, cell-cycle-confound, negative-result, falsification]
---

# The rescued "structure beyond dispersion" does not survive cell-cycle control

## Chain of evidence

1. **Structure-driving cells are proliferating cells** (scripts/18). Leave-one-cluster-out at
   osi D14: removing one Leiden cluster drops max-H1 by 0.91 (another by 0.79); all others 0.
   DE of that cluster vs rest (recovered gene symbols) = pure cell cycle — PTTG1, UBE2S, CKS1B,
   CENPW, CDKN3, MYBL2; proliferation markers 5/6 (MKI67, PCNA, CCNB1, CDK1, BIRC5); EMT 1/19,
   persister 0/14. NOT low-quality cells (median counts 10142 > 8692).

2. **Removing cycling cells reduces the excess** (stored embedding): ratio ALL=1.41 → G1-only=1.21
   (vs random-matched 1.12). Cell cycle contributes; single-run hinted it wasn't the whole story.

3. **DECISIVE — pooled embedding, bootstrapped, with proper S/G2M regression** (scripts/21):

   | preprocessing | D0 ratio | D14 ratio | D14>D0 | D14 excess(>1) |
   |---|---:|---:|---:|---:|
   | standard | 1.01 | 1.18 | 68% | 80% |
   | **S/G2M regressed** | 1.07 | **0.99** | **31%** | **33%** |

   After cell-cycle regression the D14 excess-over-dispersion is **gone (0.99)** and the drug
   increase is **gone** (D14 > D0 in only 31% of bootstraps). The standard-pooled excess (1.18)
   is itself weaker than the stored object (1.42), confirming fragility to preprocessing.

## Why the paper's cell-cycle control missed this

The manuscript's S/G2M-regression control reported that raw max-H1 was *preserved/increased*
after regression and concluded the signal is not cell cycle. But raw max-H1 is dominated by
dispersion; regression can leave raw max-H1 high (D14 ccreg realH1 = 2.39) while the
**dispersion-controlled** excess collapses to 1.0. Controlling for dispersion (Gaussian null)
AND removing cell cycle is required — and under both, no robust drug-induced structure remains.

## Honest verdict for the primary (osimertinib cell-line) system

The drug-associated max-H1 signal is explained by **transcriptional dispersion + cell-cycle
structure**. Neither "cyclic plasticity" (cells don't traverse the loop) nor "genuine
topological reorganization beyond dispersion" (does not survive cell-cycle control) is supported.

## What remains genuinely publishable

A **methodological / cautionary** contribution, not a biological discovery:
- A control framework for TDA-on-scRNA — covariance-matched Gaussian dispersion null +
  circular-coordinate negative control + bootstrap directional testing + cell-cycle regression
  in the dispersion-controlled frame.
- A case study showing an apparent topological drug-tolerance signal dissolves under these
  controls — useful to a field (persistent homology on single-cell data) that currently lacks
  standard nulls.

## Not fully closed

- PDX (GSE243562) could not be cell-cycle-scored from its processed object (HVG lacks CC genes);
  a raw rebuild is needed to test whether the strong PDX-residual excess (ratio 1.95) is also
  cell cycle. Given the primary-system result, this is unlikely to rescue the thesis but is the
  one remaining empirical gap.
