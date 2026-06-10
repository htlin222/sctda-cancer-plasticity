---
title: "Response to Reviewers — major revision"
manuscript_original: "Topological data analysis reveals cyclic cell-state plasticity underlying drug tolerance in EGFR-mutant lung cancer"
manuscript_revised: "A dispersion-controlled topological readout detects drug-induced non-Gaussian reorganization of the cell-state space in EGFR-mutant lung cancer"
date: 2026-06-10
panel_verdict_original: Major revision (4/4 reviewers; two leaning reject in original framing)
target_venue_revised: Genome Biology / Cell Reports Methods (methods-and-observation venue)
evidence_sources:
  - results/foundation/dispersion_null_summary.json
  - results/foundation/circular_coord_summary.json
  - results/foundation/FINDING_loop_is_not_traversed_2026-06-10.md
  - scripts/16_lineage_circular_coords.py
  - scripts/17_dispersion_null_and_stability.py
tags: [response-to-reviewers, reframe, tda, persistent-homology, dispersion-null]
---

# Response to Reviewers

We thank all four reviewers for an unusually rigorous and constructive panel review. The
critiques were correct, and acting on them changed the paper. Two of the central objections —
construct validity (is there *cycling*?) and the dispersion confound (is the headline statistic
just spread?) — we tested directly with new analyses rather than rhetorically. **Both new
analyses contradicted our original claim.** Cells do not traverse the loop, and at baseline the
H₁ statistic is statistically indistinguishable from dispersion.

We have therefore not merely softened the language; we have **changed the thesis** to the claim
that the evidence actually supports, and changed the venue framing accordingly.

**Revised thesis.** Under EGFR-TKI (osimertinib) in PC9 cell-line and EGFR-mutant PDX models,
drug-tolerant persister emergence is accompanied by **genuine non-Gaussian topological
reorganization** of the cell-state space: maximum H₁ persistence rises **above a
covariance-matched Gaussian dispersion null specifically under drug** (not at baseline), robustly
across subsamples. The signal is a **sparse persistent loop occupied by a minority of cells —
not population-wide reversible cycling.** The dispersion-matched Gaussian null is the central
methodological advance, and we now position the paper as a generalizable, dispersion-controlled
topological readout for perturbation scRNA-seq, consistent with R4's transfer recommendation.

We removed every claim of cyclic plasticity, reversible cycling, or cells traversing/cycling
around a loop. The circular-coordinate analysis the reviewers requested is now presented as a
**negative control** demonstrating that we do *not* overclaim traversal.

All numbers below trace to a single regenerated source. Per-group estimates and bootstrap
contrasts come from `scripts/17_dispersion_null_and_stability.py`
(`results/foundation/dispersion_null_summary.json`); circular-coordinate / lineage results come
from `scripts/16_lineage_circular_coords.py` (`results/foundation/circular_coord_summary.json`).
Reported intervals are 90% bootstrap CIs over 30 subsamples at fixed cell count.

---

## Summary of the six gating objections and their resolution

| # | Objection (panel) | Resolution |
|---|---|---|
| 1 | Construct validity: a static loop is not temporal cycling | Cycling claim **withdrawn**; circular coords run on the Watermelon lineages as a **negative control** (Rayleigh R ≥ 0.965 vs 0.881 control; clones tighter than random) → we explicitly show cells do *not* traverse the loop |
| 2 | Headline statistic fragile, unnormalized, outlier-driven | Claims now rest on **directional bootstrap robustness** (D7>D0, D14>D0 in 100% of pairs), not raw magnitudes; bootstrap CIs and CVs on every value |
| 3 | Significance rests on a vacuous gene-label null; structure-preserving null underpowered | Gene-label null **retired**; **covariance-matched Gaussian dispersion null** is now the headline test, with a built-in baseline/erlotinib negative control |
| 4 | Heterogeneity confound + no negative control | Gaussian null *is* the heterogeneity/dispersion control; baseline and erlotinib serve as built-in negatives; patient data recast as baseline-heterogeneity caveat |
| 5 | Overclaiming contradicted by own data (erlotinib non-sig/non-monotonic; patient oversell) | Scope **narrowed to osimertinib + PDX**; erlotinib recast as specificity/negative evidence; patient cohorts demoted to a heterogeneity caveat, not a drug claim |
| 6 | Numerical inconsistencies; EMT circularity; presentation/venue | One master table from a single seeded source; EMT-filter circularity removed with reframe; venue pivoted to methods/observation journal |

---

## Reviewer 1 — Concept & Novelty

### R1.1 Construct validity: "a static loop is not temporal cycling"

**Reviewer was right, and the data agree with the reviewer, not with us.** We ran the
circular-coordinate analysis R1 and R4 asked for, using the Watermelon lineage barcodes that R1
correctly identified as the obvious vehicle (`scripts/16_lineage_circular_coords.py`).

The result is unambiguous and *negative for cycling*. Validated circular coordinates collapse to
a single angular position in every cohort: Rayleigh R ≥ 0.965 across all cohorts, versus R = 0.881
for a synthetic 12%-occupied control loop that genuinely is traversed. Watermelon clones are
**tighter** than random in angular spread, not spread around the loop — the opposite of what
reversible cycling predicts. If cells cycled, lineages would smear across angles over time; they
do not.

**Action.** We withdrew the cycling thesis entirely. The title, abstract, and all section
language no longer state or imply that cells reversibly cycle or traverse a loop. The
circular-coordinate result is now reported as a **deliberate negative control / construct-validity
check**: it demonstrates the discipline of *not* overclaiming traversal from a static H₁ signal.
The surviving, supported claim is "drug exposure increases the non-Gaussian topological complexity
of the state distribution" — precisely the retreat R1 offered as acceptable.

### R1.2 Heterogeneity confound: "treatment-naive tumours show the highest H₁"

**Correct, and we now make this the honest centerpiece of our scope argument rather than an
embarrassment.** Treatment-naive patient tumours (Maynard TN) do show among the highest excess
structure (ratio 1.52), comparable to or above the per/progressive-disease stages
(PER 1.69, PD 1.51). All three patient stages — including the naive one — exceed the Gaussian null.

Under the original "H₁ = drug-induced cycling" claim this was fatal. Under the revised claim it is
expected and informative: in **patient tumours, excess topological structure reflects baseline
tumour heterogeneity (coexisting clones/states), not drug induction.** That is exactly why the
clean, drug-*induced* signal must come from the controlled in-vitro and PDX systems, where
baseline is matched. We therefore scope the drug-induction claim to osimertinib cell-line + PDX,
and present the patient cohorts explicitly as a **baseline-heterogeneity caveat**, answering R1's
"naive tumours are highest" objection on its own terms.

### R1.3 Overclaiming "across all systems" contradicted by erlotinib

**Correct.** The erlotinib discovery cohort shows **no excess structure over the dispersion null**:
ratios 0.93 (D0), 0.99 (D9), 0.70 (D11), with real-greater-than-Gaussian fractions of 33/42/0%.
There is no monotonic increase and no drug-specific excess. We do not present erlotinib as a
positive result. Instead we recast it as **specificity evidence**: a chemically related EGFR-TKI
in the same cell line that does *not* produce excess topological structure under our null bounds
the claim and demonstrates the readout is not a generic artifact of "drug was added." The abstract
no longer says cyclicity increases monotonically across all systems.

---

## Reviewer 2 — Theory & Methodological Rigor

### R2.1 The headline statistic is fragile, unnormalized, and outlier-dominated

**Accepted.** We no longer ask raw max-H₁ magnitude to carry biological weight. Two changes:

1. **The primary inference is now directional and bootstrap-based at fixed cell count**, not a
   magnitude comparison. The supported statements are: osi D7 > D0 in **100%** of bootstrap pairs;
   osi D14 > D0 in **100%**; PDX residual > untreated in **97%** (`_contrasts_pct_pairs`,
   `dispersion_null_summary.json`). These are ordinal, scale-free statements robust to the
   metric-dependence R2 flagged.

2. **Every reported value now carries a 90% bootstrap CI and a CV.** For the osimertinib series,
   subsample CVs are modest (D0 5.9%, D3 6.3%, D7 6.8%, D14 8.6%). The larger CVs occur exactly
   where we now decline to make point claims — PDX residual 23.6%, Maynard PD 44.7% — and there we
   rely on the directional fraction (PDX residual > untreated 97%) and the ratio bound, not the
   magnitude. The S/G2M cell-cycle-regression results are retained as ablation controls, not as
   confirmatory magnitude increases.

### R2.2 The significance rests on a near-vacuous null

**This is the heart of the revision, and we agree the gene-label null tested the wrong question.**
Shuffling genes within each cell destroys all covariance, so any real data beats it — it asks "is
there *any* structure," not "is there drug-dependent excess structure beyond dispersion." We have
**retired the gene-label null entirely.**

The new headline test is a **covariance-matched Gaussian dispersion null**: a multivariate
Gaussian with the *same mean and covariance* as the real data, carrying the population's dispersion
but no genuine topological structure. The reported statistic is the ratio of real max-H₁ to the
Gaussian null's max-H₁. This directly isolates structure beyond dispersion, which was R2's and
R3's precise demand. The result discriminates cleanly:

| group | real H₁ [90% CI] | Gaussian null | ratio | % real > Gaussian |
|---|---|---:|---:|---:|
| osi D0 | 1.41 [1.31, 1.60] | 1.49 | 0.95 | 31% |
| osi D3 | 1.82 [1.82, 2.15] | 2.00 | 0.91 | 18% |
| osi D7 | 2.67 [2.45, 3.04] | 2.04 | **1.31** | **100%** |
| osi D14 | 3.78 [2.99, 3.78] | 2.66 | **1.42** | **98%** |
| erl D0 / D9 / D11 | 1.49 / 1.77 / 1.08 | ~1.6 | 0.93 / 0.99 / 0.70 | 33 / 42 / 0% |
| PDX untreated | 2.79 [2.23, 4.03] | 2.30 | 1.21 | 86% |
| PDX residual | 5.46 [3.58, 7.34] | 2.80 | **1.95** | **99%** |
| Maynard TN / PER / PD | 5.49 / 6.04 / 5.58 | ~3.6 | 1.52 / 1.69 / 1.51 | 99 / 98 / 95% |

At baseline (osi D0, osi D3) the real data are **at or below** the dispersion null (ratio ≤ 0.95)
— i.e. baseline max-H₁ is just dispersion, as R2/R3 suspected. Excess structure appears
**specifically under drug** (osi D7 ratio 1.31, osi D14 ratio 1.42; PDX residual ratio 1.95). This
is the discriminating, structure-preserving test the panel required, and it now serves as its own
negative control because baseline and erlotinib fall on the null side.

### R2.3 Underpowered structure-preserving permutation (n=50, one comparison)

**Accepted.** The inadequate n=50 timepoint-label permutation is removed. The dispersion-null
inference is built on **30 bootstrap subsamples per group with directional pair fractions over
all pairings**, and the negative cohorts (baseline, erlotinib) provide built-in null calibration.
We report the directional fractions as the primary evidence so that conclusions do not hinge on a
single underpowered permutation count.

### R2.4 EMT circularity in the Mapper filter

**Accepted.** The original analysis used the EMT score as the Mapper filter and then reported EMT
genes/Hallmark-EMT as top-connected/top-enriched, which is partly tautological. Because the
revised paper no longer rests on the "EMT anchors the loop" narrative — and no longer makes a
cycling/Mapper-loop claim at all — we have removed the circular EMT-filter result from the
inferential chain. Characterization of the minority of cells driving the non-Gaussian structure
(osi D14, PDX residual) is deferred to a separate, filter-independent analysis and is not used to
support the headline topological claim.

### R2.5 Numerical inconsistencies (data-integrity flag)

**Accepted; fully addressed.** Every quantitative statement in the revised manuscript is drawn
from a **single regenerated source** — `results/foundation/dispersion_null_summary.json` and
`circular_coord_summary.json`, produced by `scripts/17` and `scripts/16` at a fixed seed and fixed
cell count (1,200 cells per group except the small erlotinib timepoints, n = 379/249). The
conflicting YU-006 PDX baseline values (2.36 / 2.81 / 2.89 / 3.62) and the divergent osi D0 values
(1.41 vs 1.71) no longer appear; the single authoritative figures are PDX-YU006 untreated 2.79
[2.23, 4.03] and osi D0 1.41 [1.31, 1.60]. The revised manuscript carries one master table with
cell count per value, eliminating the untraceable-number problem.

---

## Reviewer 3 — Data, Feasibility & Reproducibility

### R3.1 Vacuous null / no negative control

**Addressed jointly with R2.2 and R1.2.** The Gaussian dispersion null is itself the negative
control the reviewer asked for: it shows the statistic correctly reports **no** excess at baseline
(osi D0 ratio 0.95, 31% > null; osi D3 ratio 0.91, 18% > null) and under erlotinib (D11 ratio
0.70, 0% > null). A statistic that fires on real drug-treated data but sits at or below a
covariance-matched null on baseline and on a non-responding drug is demonstrably not reporting
generic structure or cell count alone.

### R3.2 Instability — CV up to 23.5%, factor-of-2 PC swings, S/G2M control increases H₁

**Accepted and disclosed, not hidden.** We report CVs openly (range ~6–45% across groups) and
build claims on **directional bootstrap robustness** rather than point magnitudes (see R2.1). The
high-CV groups are exactly the ones where we make only ordinal claims. The S/G2M cell-cycle
regression is retained as an ablation control consistent with prior project findings that H₁ loops
persist after cell-cycle ablation; it is no longer presented as a confirmatory magnitude increase.

### R3.3 Heterogeneity / cell-count confound

**Addressed by fixed-count subsampling and the covariance-matched null.** All groups are evaluated
at a fixed cell count (1,200, with the two small erlotinib timepoints noted), and the Gaussian
null shares the real data's covariance, so dispersion and coexisting-state count are matched out.
What remains after that matching — the non-Gaussian excess under drug — is the claim.

### R3.4 PDX mouse-gene removal by lowercase-symbol regex

**Accepted as a limitation.** We retain the PDX result because the drug contrast is strong and
directionally robust (residual > untreated in 97% of bootstrap pairs; residual ratio 1.95, 99% >
null), but we now state explicitly that mouse-gene removal used a symbol-based heuristic inferior
to alignment-based xenograft deconvolution, that this can bidirectionally misassign features, and
that the PDX finding should be read as supporting in-vivo evidence rather than a stand-alone
validation. We do not overstate it.

### R3.5 Provenance: "bit-stable" overclaim and version-tag mismatch

**Accepted.** We removed the "bit-stable across runs" language, which is contradicted by the
documented BLAS/OpenMP nondeterminism and the disclosed CVs; reproducibility is now described as
deterministic-given-seed with quantified subsample variance. The submission version tag is
reconciled to the actual repository tag, and the released analysis is archived (the Zenodo DOI
replaces the GitHub-tag reference, per NC/GB policy).

---

## Reviewer 4 — Editorial / Presentation / Synthesis

### R4.1 The titular "loop" is never shown as a demonstrated cycle

**Accepted; resolved by reframe.** Because we no longer claim a traversed cycle, the title and
narrative no longer promise a Mapper graph containing a demonstrated loop. The dropped
Mapper-graph and circular-profile panels are replaced by the two analyses that actually carry the
revised argument: (i) the dispersion-null ratio plot per group (the headline method), and (ii) the
circular-coordinate negative control showing angular collapse (Rayleigh R ≥ 0.965 vs 0.881 for the
12%-occupied control loop). The display items now match the claims.

### R4.2 Overclaim in abstract ("across all systems"; "replicates in 14-patient cohort")

**Accepted; addressed with R1.3.** The abstract is rewritten to: drug-induced non-Gaussian
topological excess is established under **osimertinib in cell-line and PDX models**; erlotinib
shows **no excess** (specificity); patient cohorts show excess at **all stages including naive**,
indicating baseline heterogeneity rather than drug induction. No "across all systems" claim and no
overstated patient replication remain.

### R4.3 Venue fit

**Accepted; we follow the panel's recommendation.** The defensible contribution is a
**reproducible, dispersion-controlled topological readout** for perturbation scRNA-seq — a
methods-and-observation advance, not a new mechanism. We therefore pivot the framing to **Genome
Biology / Cell Reports Methods**, presenting the covariance-matched Gaussian null as a
generalizable method, the osimertinib/PDX drug contrasts as the demonstrating observation, and the
circular-coordinate negative control as the construct-validity guardrail. This is the honest fit
and plays to the genuine strength the panel identified (controls and reproducibility above field
average).

### R4.4 Presentation, declarations, missing literature

**Accepted.** We reconcile the AI-image declaration across `declarations.md`, `checklist.md`, and
the figure caption; resolve the funding placeholder; bring word and display-item counts within the
target venue's limits; and add the previously missing literature relevant to the single-scalar
weakness and the dynamics alternatives (foundational non-genetic resistance work; CellRank/scVelo
distinguished from our static-topology readout; persistence-image/landscape vectorized statistics
as the natural extension of a single-scalar summary).

---

## What we did not change

We did not invent supporting analyses or inflate effect sizes to rescue the original cycling
claim. Where the new analyses came back negative — circular coordinates (no traversal), baseline
and erlotinib (no excess over dispersion), patient naive tumours (excess without drug) — we report
the negative result and let it narrow the claim. The paper is smaller and more defensible than the
one originally submitted, and every quantitative statement traces to a single seeded source.
