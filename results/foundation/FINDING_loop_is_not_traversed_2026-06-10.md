---
title: "Foundation finding — the H1 loop is geometrically concentrated, not traversed by the cell population"
date: 2026-06-10
analysis: scripts/16_lineage_circular_coords.py
dataset: GSE150949 PC9 osimertinib (Watermelon lineage), processed h5ad, top-30 PCA
method: dreimac CircularCoords (DSPVJ), prime=47, n_landmarks=400
status: first foundation analysis; robust to perc; needs extension to other cohorts/PC counts
bearing_on: peer-review construct-validity objection (reviewer #1/#3/#4)
tags: [foundation, construct-validity, circular-coordinates, tda, falsification-test]
---

# Foundation finding: the H₁ loop is real topology, but the population does **not** traverse it

## What was tested

The manuscript's central interpretation is that drug-tolerant cells **reversibly cycle**
among states, detected as rising max-H₁ persistence. max-H₁ is computed on a **static
snapshot**. This analysis asks the construct-validity question directly: *if* there is a
genuine cell-state cycle, cells (and clonal lineages) should be distributed **around** the
H₁ loop, not piled at one spot.

We assigned each cell a circular coordinate θ ∈ [0,1) from the most-persistent H₁ cohomology
class using dreimac (validated DSPVJ implementation), on the same top-30 PCA used in the paper.

## Method is validated (two checks)

1. **Reproduces the paper.** dreimac's H₁ persistence at D14 = **3.776**, matching the
   manuscript's osimertinib D14 = 3.78. Same topological feature.
2. **Positive controls pass.**
   - Clean noisy circle → **R = 0.124, coverage = 1.00** (cells fill the circle, as they should).
   - Dense blob with only **12 %** of cells on a thin loop → **R = 0.881, coverage = 0.83**.

   (R = mean-resultant length: R→0 = spread around circle, R→1 = all at one angle.
   coverage = fraction of 12 circular sectors occupied.)

## Result (osimertinib, all four timepoints)

| Timepoint | H₁ persistence | circle coverage | Rayleigh R |
|-----------|---------------:|----------------:|-----------:|
| D0  | 1.57 | 0.25 | 0.998 |
| D3  | 2.22 | 0.33 | 0.999 |
| D7  | 2.45 | 0.50 | 0.996 |
| D14 | 3.78 | 0.33 | 0.997 |

- **R ≈ 0.997–0.999 at every timepoint** — *more concentrated than a cloud with only 12 % of
  cells on its loop.* The population sits at essentially one angle; it does not occupy the circle.
- **Coverage 0.25–0.50** — only 3–6 of 12 sectors filled.
- **Robust to dreimac `perc` (0.3→0.9):** R stays 0.997–0.998, coverage unchanged. Not a
  parameter artifact.

## Clonal lineage test (Watermelon barcodes)

Within each timepoint, for clones with ≥3 co-sampled cells, we measured angular circular
variance (CV: 0 = clone at one angle, 1 = clone spread around loop) vs. size-matched random
cell groups (1000 permutations).

- Real clones have **CV ≈ 0** at every timepoint — clonemates sit at one angle.
- At **D14**, real clones are **more concentrated than random** (p = 0.036 for "more spread"
  is *rejected*; clones are tighter than random, p = 0.965). This is the **opposite** of the
  "clones explore a state cycle" prediction.

## Cross-cohort confirmation (all five systems, validated metric, 1500-cell subsamples)

Repeating the circular-coordinate spread on each cohort's headline group:

| Cohort group | R (1=one angle) | coverage | % cells off dominant wedge |
|--------------|----------------:|---------:|---------------------------:|
| erlotinib D9 (paper's "discovery") | 0.992 | 0.33 | 1.3 % |
| osimertinib D14 | 0.997 | 0.33 | 0.6 % |
| PDX YU-006 residual | 0.965 | 0.83 | 6.1 % |
| Kim LUAD naive | 0.987 | 0.67 | 48 %* |
| Maynard PD | 1.000 | 0.25 | 0.1 % |

Reference: the positive control with **12 %** of cells on a thin loop gave R = 0.881.
**Every cohort has R ≥ 0.965 > 0.881** — i.e. in all five systems the most-persistent loop is
occupied by *fewer* cells than that 12 % control, and the population does not circumnavigate it.
*Kim's 48 % off-dominant-wedge with R = 0.987 is an *arc* across adjacent sectors (a partial
spread in one region), not a traversed circle.

## Outlier-robustness (important correction)

A separate check trimmed the top-1 % most-peripheral cells (by 15-NN distance) and recomputed
max-H₁. It did **not** collapse in 7/9 groups (only osi-D14 −21 %, Maynard-TN −17 %). So the
statistic is **not** a pure "single extreme-outlier" artifact, and an earlier note in this file
inferring "the loop = 4 cells" from the *cohomology cocycle support* was an over-read — a
cocycle's support is the dual cut, not the count of cells on the loop. **Retracted.** The
robust, validated claim is the circular-coordinate one above: the loop is occupied by a small
minority and is not traversed, but it is a reproducible sparse feature, not one freak cell.

## Interpretation

The persistent H₁ bar the paper reports is a **real, reproducible topological feature**, but
across all five cohorts it is **occupied by a small minority of cells (≲10 %, mostly ≲2 %) and
is not circumnavigated** — the circular coordinate collapses to one angle/arc, and in
osimertinib clonal lineages are *tighter* than random, not spread around the loop. Neither the
population nor the clones traverse a cycle in state space.

**This corroborates the peer-review construct-validity objection and decisively answers the
"can we save the cycling interpretation by finding loop-resident cells?" question: no.** In no
cohort does a substantial cell population reside on or traverse the loop, so there is no
loop-resident program to anchor a "cycling" claim. A rising max-H₁ does not demonstrate
"reversible inter-state cycling"; it tracks a sparse minority feature.

## What this does NOT claim (caveats / next checks)

- It does **not** dispute that max-H₁ rises with drug exposure (it does; reproduced here:
  osi D14 = 3.78, etc.).
- It does **not** claim the loop is one freak outlier (trim test; see correction above).
- Each cell is sequenced once; "traversal" is tested via the clonal-spread proxy, which is
  the strongest available but indirect.
- Top-30 PCA only; a PC-count sweep would further harden the cross-cohort claim.
- A concentrated cloud with a minor persistent arc is consistent with "one dominant state +
  a small sub-population arc," i.e. max-H₁ may track a minority feature rather than
  population-wide plasticity.

## Bearing on the manuscript

Supports reframing from **"cyclic / reversible plasticity"** toward a defensible claim about
**a localized topological feature / increased complexity** under drug exposure — unless a
follow-up shows the loop is genuinely cell-traversed (e.g. loop-resident cells form a
coherent expression program that varies smoothly around θ). This is the first analysis that
touched the scientific foundation rather than the manuscript's packaging.
