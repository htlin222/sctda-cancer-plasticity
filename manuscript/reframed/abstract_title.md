# Reframed title and abstract

Scope: osimertinib-treated PC9 cell line and EGFR-mutant PDX. Erlotinib = specificity
(no-excess) cohort; patient tumours = baseline-heterogeneity caveat, not a drug claim.
All quantitative statements trace to the locked evidence table; no new numbers.

## Title options

1. **A dispersion-controlled topological signature of drug-tolerant cell-state reorganization in EGFR-mutant lung cancer**

2. **Drug tolerance under EGFR inhibition emerges as non-Gaussian topological reorganization of the single-cell state space**

3. **A covariance-matched null reveals drug-specific topological structure in the cell-state space of EGFR-mutant persister cells**

Preferred: Option 1 (names the methodological advance — dispersion control — and the
biological observation — drug-tolerant reorganization — without any cycling/traversal claim).

## Abstract (~150 words)

Drug-tolerant persister cells survive EGFR-targeted therapy and seed relapse in lung
cancer, but whether their emergence reshapes the global geometry of cell-state space has
been difficult to test rigorously. A central obstacle is that topological summaries of
single-cell data scale with mere dispersion, inflating apparent structure. Here we
introduce a covariance-matched Gaussian null that isolates genuine non-Gaussian structure
from dispersion, and apply it to longitudinal single-cell RNA-seq of EGFR-mutant lung
cancer. In osimertinib-treated PC9 cells, the maximum H1 persistence rises above the
dispersion-matched null specifically under drug (ratio 1.31 at day 7, 1.42 at day 14;
day-7 and day-14 exceed baseline in 100% of bootstrap pairs), with no excess at baseline.
The same drug-specific excess appears in EGFR-mutant patient-derived xenografts (residual
ratio 1.95). The structure is a sparse persistent loop occupied by a minority of cells;
circular-coordinate analysis confirms cells do not traverse it. Erlotinib shows no excess,
underscoring specificity.

## Why no cycling claim survives

The original abstract claimed "reversible cyclic plasticity" and that "cyclic cell-state
structure increases monotonically with drug exposure across all systems." Both are removed:
- Cycling/traversal is recast as a built-in negative control. The persistent loop is sparse
  and minority-occupied; circular coordinates collapse (Rayleigh R >= 0.965), so cells do
  not traverse it. The abstract states this explicitly.
- "Monotonic across all systems" is dropped. Scope is narrowed to osimertinib cell line and
  PDX; erlotinib is presented as a no-excess specificity cohort; patient tumours (excess at
  all stages including naive) are not claimed as drug induction.
- The headline is the methodological advance: a covariance-matched Gaussian null that
  separates real topological structure from dispersion.

## Number provenance (locked)

- osi D7 ratio 1.31 (real max-H1 2.67 [2.45,3.04] vs gauss 2.04), 100% bootstrap real>gauss
- osi D14 ratio 1.42 (real 3.78 [2.99,3.78] vs gauss 2.66), 98% bootstrap real>gauss
- osi D0 ratio 0.95 (real 1.41 [1.31,1.60] vs gauss 1.49) — baseline = dispersion
- Directional: osi D7>D0 and D14>D0 in 100% of bootstrap pairs
- PDX residual ratio 1.95 (real 5.46 [3.58,7.34] vs gauss 2.80), 99%; >untreated in 97%
- Erlotinib D0/D9/D11 ratios 0.93/0.99/0.70 — no excess
- Circular coordinates: Rayleigh R >= 0.965 (vs 0.881 for 12%-occupied control loop) — no traversal
