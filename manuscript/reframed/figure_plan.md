# Reframed figure plan — dispersion-controlled topological reorganization under EGFR-TKI

Target venue: methods + observation journal (Genome Biology / Cell Reports Methods).
Style: 4 main figures + 1 optional, 180 mm max width, 7–8 pt sans, 300 dpi vector PDF.

**Thesis carried by the figures (every panel must serve this):**
Under EGFR-TKI (osimertinib) in PC9 cell line and EGFR-mutant PDX, drug-tolerant persister
emergence is accompanied by **genuine non-Gaussian topological reorganization** of cell-state
space: maximum $H_1$ persistence rises **above a covariance-matched Gaussian dispersion null
specifically under drug, not at baseline**, robustly across bootstrap subsamples. The signal is
a **sparse persistent loop occupied by a minority of cells**, explicitly **not** population-wide
cycling — a negative control proves cells do not traverse the loop. The dispersion-matched
Gaussian null is the central methodological advance.

**Hard constraints applied throughout (see HARD RULES in the rescue brief):**
- No "cyclic plasticity", "reversible cycling", "cells traverse/cycle a loop", no donut/cycle
  concept figure. The circular-coordinate panel is framed as a NEGATIVE control against the
  traversal overclaim.
- Scope = osimertinib cell line + PDX. Erlotinib = no-excess / specificity evidence. Patients =
  baseline-heterogeneity caveat, not a drug claim.
- Every number below traces to the locked evidence table. No new numbers.

---

## Figure 1 — The dispersion-matched null and the max-$H_1$ readout (METHOD)

**Message:** "Maximum $H_1$ persistence alone conflates real structure with dispersion; the
covariance-matched Gaussian null separates the two. This is the method."

- **(a) Schematic of the readout.** A point cloud in cell-state (PCA) space; the most-persistent
  $H_1$ generator is highlighted; a Vietoris–Rips filtration cartoon shows the loop being born
  and dying, defining *maximum $H_1$ persistence*. No "cycle/traversal" arrows — purely the
  filtration-scale definition.
- **(b) The dispersion confound, illustrated.** Two synthetic clouds with **identical mean and
  covariance**: one structureless Gaussian, one with an embedded loop. They yield similar
  max-$H_1$, motivating why a covariance-matched null — not a gene-label shuffle — is required.
- **(c) Null construction.** Diagram: from each real subsample, draw a multivariate Gaussian
  with the **same mean and covariance**, recompute max-$H_1$; repeat to build the null. Report
  the per-cohort ratio = real / Gaussian-null and the bootstrap fraction (real > null).
- **(d) Baseline = dispersion (the null behaves as designed).** osi D0 real max-$H_1$
  1.41 [1.31, 1.60] vs Gaussian 1.49, **ratio 0.95**, real > null in only **31%** of bootstraps;
  osi D3 1.82 [1.82, 2.15] vs 2.00, **ratio 0.91**, **18%**. At baseline the real data does not
  exceed its own dispersion null.

*Panel-d numbers: osi D0/D3 from dispersion_null_summary.json.*

---

## Figure 2 — Drug-specific emergence of excess structure in PC9 (HERO)

**Message:** "Under osimertinib, and only under drug, real max-$H_1$ rises above the
dispersion-matched null — robustly across subsamples."

- **(a) HERO panel — real vs Gaussian-null across the osimertinib time course.** Paired
  bars/dumbbells per timepoint, real max-$H_1$ with 90% bootstrap CI against the Gaussian-null
  median:
  - D0 1.41 [1.31, 1.60] vs 1.49 — ratio 0.95 (31% > null)
  - D3 1.82 [1.82, 2.15] vs 2.00 — ratio 0.91 (18% > null)
  - D7 2.67 [2.45, 3.04] vs 2.04 — **ratio 1.31 (100% > null)**
  - D14 3.78 [2.99, 3.78] vs 2.66 — **ratio 1.42 (98% > null)**
  Annotate the crossover between D3 and D7 where real first exceeds the null. This panel is the
  paper's central evidence.
- **(b) Excess-structure ratio vs day.** Line of ratio (real/null) across D0→D14, with a
  reference line at ratio = 1; the curve crosses 1 between D3 and D7 and rises to 1.42 at D14.
- **(c) Directional bootstrap robustness.** Bar/text panel: osi **D7 > D0 in 100%** and
  **D14 > D0 in 100%** of paired bootstrap subsamples. Claims rest on this directional
  robustness, not on unstable point estimates (per-cohort CV 6–9% disclosed in caption).

*All numbers from dispersion_null_summary.json (osi D0/D3/D7/D14, _contrasts_pct_pairs).*

---

## Figure 3 — The loop is sparse and not traversed (NEGATIVE CONTROL / construct validity)

**Message:** "The excess structure is a sparse loop occupied by a minority of cells. Cells do
NOT circumnavigate it — we explicitly do not claim cycling."

- **(a) Circular-coordinate collapse.** Validated circular coordinates (dreimac DSPVJ,
  prime = 47, 400 landmarks) on the osi D14 / PDX-residual most-persistent loop. Angle histogram
  (or rose plot) showing cells concentrated at essentially one angle. Rayleigh **R = 0.997**
  (osi D14), and **R ≥ 0.965 in all five cohorts**.
- **(b) Negative-control benchmark.** Side-by-side with a synthetic positive control where
  **12% of cells lie on a thin loop**, giving **R = 0.881**. Every real cohort has **R ≥ 0.965 >
  0.881** — i.e. occupied by *fewer* cells than the 12% control. Interpretation banner: "the
  population sits at one angle; it does not occupy the circle → no traversal."
- **(c) Clonal-lineage test (Watermelon barcodes).** Clones are *tighter* than the random null
  (concentration not rejected; spread rejected, p = 0.965) — the **opposite** of the cycling
  prediction. Cells of a lineage do not spread around the loop.
- **(d) Sparse-occupancy summary.** Small table/strip: across cohorts the loop is occupied by a
  small minority (≲10%, mostly ≲2%) of cells. Frames the surviving claim: a *sparse* persistent
  loop, not population-wide cycling.

*Numbers from circular_coord_summary.json (R = 0.997 osi D14) and
FINDING_loop_is_not_traversed_2026-06-10.md (control R = 0.881 at 12%; R ≥ 0.965 all cohorts;
Watermelon p = 0.965).*

---

## Figure 4 — In-vivo replication in PDX and honest scope (REPLICATION + SPECIFICITY)

**Message:** "The drug-induced excess structure replicates in EGFR-mutant PDX in vivo;
erlotinib shows no excess (specificity), and patient tumors carry baseline heterogeneity (not a
drug claim)."

- **(a) PDX in-vivo replication (hero of this figure).** Real vs Gaussian-null paired bars,
  90% CI:
  - YU-006 untreated 2.79 [2.23, 4.03] vs 2.30 — ratio 1.21 (86% > null)
  - YU-006 residual 5.46 [3.58, 7.34] vs 2.80 — **ratio 1.95 (99% > null)**
  Annotate **PDX residual > untreated in 97%** of paired bootstraps. Independent in-vivo system
  reproduces the drug-on excess.
- **(b) Erlotinib = no excess (specificity).** Real vs Gaussian-null for erl D0/D9/D11:
  1.49 / 1.77 / 1.08 vs ~1.6/1.80/1.55, ratios **0.93 / 0.99 / 0.70**, never exceeding the null
  (D11 0% > null). Framed as a specificity control: the readout does not fire for a cohort
  without the osimertinib drug-on signal, strengthening that the D7/D14/PDX-residual effect is
  not a generic artifact.
- **(c) Patient baseline-heterogeneity caveat.** Maynard TN / PER / PD real vs null:
  5.49 / 6.04 / 5.58 vs ~3.6, ratios **1.52 / 1.69 / 1.51** — excess at **all** stages including
  **treatment-naive (TN)**. Explicit caption: in patient tumors, excess topological structure is
  present at baseline and reflects **tumor heterogeneity, not drug induction**; therefore the
  clean drug-induction claim is confined to the controlled cell-line and PDX systems.

*Numbers from dispersion_null_summary.json (PDX-YU006 untreated/residual, erl D0/D9/D11, Maynard
TN/PER/PD) and _contrasts_pct_pairs (PDX residual > untreated 97%).*

---

## Figure 5 (OPTIONAL / supplement-promotable) — Cell-cycle robustness of the drug-on signal

**Message:** "The osimertinib D7/D14 and PDX-residual excess survives removal of the dominant
cell-cycle variance axis."

- **(a)** Tirosh cell-cycle gene ablation: max-$H_1$ before/after with persistence ratio
  ≥ 0.71 across timepoints (existing control; loops are not cell-cycle artifacts).
- **(b)** Stringent S/G2M regression on 1,000-cell subsamples: PDX YU-006 untreated 2.89→3.67
  (ratio 1.27) and residual 4.65→5.96 (ratio 1.28); the residual > untreated ordering is
  preserved after removing the cell-cycle axis.

*Uses existing cell-cycle control numbers already in tab:cc_controls. Promote to supplement if
the 4-figure layout is preferred for a Methods venue.*

---

## Mapping from the old (retracted) plan

| Old main figure | Fate |
|---|---|
| Fig 1 study design + velocity/PAGA "tree" | Folded into Fig 1a schematic + Methods; UMAP/velocity demoted to supplement |
| Fig 2 "emerging loops" (max-$H_1$ + gene-label permutation null) | Replaced by Fig 1 (Gaussian null) + Fig 2 (hero); gene-label null retired |
| **Concept figure: donut / reversible cycling** | **DROPPED entirely (overclaim)** |
| Fig 3 Mapper "loop is real" | Dropped as a positive structure claim; Mapper at most a supplement |
| Fig 4 cell-cycle ablation | Demoted to optional Fig 5 / supplement |
| Fig 5 "reversible EMT/MET cyclic transition", polar gene profiles, cyclic-state model | **DROPPED** (cycling claim + circular gene polar plots imply traversal) |
| Fig 6 + master "moves toward patient topology" | Dropped; patient data recast as baseline-heterogeneity caveat in Fig 4c |

---

## Cross-figure conventions

- **Single source of truth:** every max-$H_1$, CI, ratio, and bootstrap % comes from
  `results/foundation/dispersion_null_summary.json` and `circular_coord_summary.json`
  (one seed/source) to eliminate the numerical-inconsistency objection.
- **Always pair real with its Gaussian null** wherever max-$H_1$ appears; never show a bare
  max-$H_1$ value as evidence of structure.
- **Always show 90% bootstrap CIs** on real values and report ratio + % (real > null).
- **Color:** real = dark teal; Gaussian null = neutral grey; ratio > 1 highlighted, ratio ≤ 1
  muted. Drug-on conditions warm-accented only when ratio > 1 with ≥95% bootstrap support
  (osi D7, osi D14, PDX residual).
- **No cyclic/traversal language or iconography in any caption, axis label, or legend.**
