# A control framework for topological data analysis of perturbation single-cell RNA-seq: dispersion-matched nulls, traversal tests, and a drug-tolerance case study

**Short title:** Structure-aware controls for topological analysis of scRNA-seq

**Authors:** H.-T. Lin, Y.-H. Tu

**Target journals:** Bioinformatics · GigaScience · Genome Biology (Method) · Cell Reports Methods

---

## Abstract

Persistent homology is increasingly applied to single-cell RNA-seq (scRNA-seq) to detect
"cyclic" or non-tree-like cell-state structure, with the maximum degree-1 persistence (max-$H_1$)
used as a scalar readout of cell-state plasticity. We show that this statistic, computed without
structure-aware controls, conflates genuine topology with two mundane confounds — transcriptional
dispersion and the cell cycle — and we provide a reusable four-part control framework that
separates them: (i) a **covariance-matched Gaussian null** that distinguishes real structure from
dispersion; (ii) a **circular-coordinate traversal test** that asks whether a detected "loop" is
actually occupied before any cyclic interpretation; (iii) **cell-cycle regression evaluated
within the dispersion-controlled frame**, which we show is necessary because the published-style
control (raw max-$H_1$ preserved after regression) is inadequate; and (iv) **bootstrap
directional testing**, because max-$H_1$ point estimates are subsample-unstable. On synthetic
data with known ground truth the framework is sensitive (it passes a genuine traversed loop) and
specific (it rejects benign dispersion and non-Gaussian non-cyclic structure). Applied as a
worked case study to five longitudinal EGFR-mutant lung-cancer systems, the framework decomposes
an apparently compelling, cross-system, drug-associated max-$H_1$ signal: cells never traverse the
loop (Rayleigh $R \ge 0.965$ vs $0.881$ for a 12%-occupied control), baseline max-$H_1$ equals a
dispersion-matched null in every system, the in-vitro drug excess is contributed by proliferating
cells and vanishes after cell-cycle regression, and the in-vivo (PDX) excess is genuine
non-cell-cycle structure that nonetheless still fails the traversal test. In no system does
max-$H_1$ evidence a cell-traversed cycle. We recommend these controls as a standard for
topological analyses of perturbation single-cell data and release them as open-source code.

## Contributions

1. A **covariance-matched Gaussian null** for max-$H_1$ that separates genuine topological
   structure from transcriptional dispersion, superseding the field-standard within-cell
   gene-shuffle null (which destroys all covariance and is beaten by essentially any real data).
2. A **circular-coordinate traversal test** establishing occupancy of a persistent loop as a
   precondition for any "cyclic" claim.
3. The observation — and correction — that the **common cell-cycle control (raw max-$H_1$
   preserved after S/G2M regression) is inadequate**, because raw max-$H_1$ is dispersion-
   dominated; the cell-cycle contribution is only visible when evaluated as a ratio to the
   dispersion null.
4. A **synthetic benchmark** characterising the framework's sensitivity and specificity, and a
   five-system cancer **case study** demonstrating a confounded signal that is reproducible,
   cross-scale, and wrong.

---

## Introduction

Topological data analysis (TDA), and persistent homology in particular, is an increasingly
popular lens on single-cell RNA-seq. Its appeal is the ability to detect closed, non-tree-like
structure (loops; degree-1 homology, $H_1$) that trajectory-inference tools, designed for
one-directional progressions, are not built to quantify. In cancer drug tolerance this is
especially attractive: lineage-tracing studies show drug-tolerant persister cells reversibly
transition among transcriptional states, a behaviour naturally described as "cyclic." A common
move is to summarise the topology with a single scalar, the maximum $H_1$ persistence
(max-$H_1$), and to read its increase under drug as a measure of cyclic plasticity.

This inherits a known hazard. Persistent homology on a noisy high-dimensional point cloud reports
$H_1$ features whether or not the population occupies or traverses a loop, and max-$H_1$ as a
scalar conflates genuine topological structure with two confounds: overall transcriptional
dispersion (the variance/spread of the embedding, which directly inflates Vietoris–Rips
persistence) and the cell cycle, which is itself intrinsically a loop (G1$\to$S$\to$G2M$\to$G1).
Cautionary methods work has repeatedly shown that single-cell analyses require structure-aware
nulls — for doublets, ambient RNA, batch effects, and over-interpreted pseudotime — and
topological summaries are no exception. Yet the nulls in common use for scTDA (within-cell
gene shuffles) do not control for dispersion, and cell-cycle controls are typically evaluated on
raw persistence rather than relative to a dispersion baseline.

Here we provide a control framework that makes both confounds explicit, validate it on synthetic
data with known ground truth, and apply it to five EGFR-mutant lung-cancer systems as a worked
case study. Our aim is not a new biological discovery but to prevent a class of false ones.

---

## Results

### The control framework

The framework takes a cell-by-PC embedding for each condition and reports three quantities plus a
decision rule. (1) **Dispersion test:** the ratio of observed max-$H_1$ to that of a multivariate
Gaussian matched to the embedding's mean and covariance (median over draws); a ratio near 1 means
the statistic is explained by dispersion, a ratio robustly $>1$ indicates structure beyond
dispersion. (2) **Traversal test:** the Rayleigh concentration $R$ of circular coordinates derived
from the most-persistent $H_1$ class; $R$ near 1 means the population sits at one angle (loop not
occupied), low $R$ means cells are distributed around it. (3) **Cell-cycle test:** the dispersion
ratio recomputed after S/G2M-score regression; loss of the excess attributes it to the cell cycle.
All three are wrapped in bootstrap subsampling, and claims rest on the fraction of subsamples in
which an ordering holds, not on point estimates. Decision rule: *structure beyond dispersion* if
the Gaussian ratio is robustly $>1$; *a traversed loop* only if additionally $R<0.6$; *cell-cycle
driven* if the excess does not survive S/G2M regression.

### The framework is sensitive and specific on synthetic ground truth

We validated the framework on point clouds where the answer is known (n=1,200 each, structure in
two dimensions embedded in 30-D with isotropic noise; **Fig. 4**). A pure Gaussian blob (no
structure) gives ratio 0.93 and $R=0.99$ — correctly **no structure** (true negative). A
uniformly traversed loop gives ratio 8.35 and $R=0.42$ — correctly a **traversed loop** (true
positive; sensitivity). A sparse non-closed arc on a blob gives ratio 1.54 but $R=0.96$ —
correctly **structure that is not traversed**, reproducing the real-data signature below. Two
well-separated Gaussian clusters — non-Gaussian but non-cyclic — give ratio 0.93, correctly **no
structure** (specificity): the Gaussian null does **not** false-positive on benign
non-Gaussianity, the key objection to a covariance-matched null. A second, non-parametric
dispersion-preserving null (independent per-PC permutation) agrees on the blob, loop, and clusters
and is more conservative on the sparse arc (ratio 0.96), confirming that the dispersion conclusion
does not rest on Gaussianity alone; we use the Gaussian null as primary (sensitive and specific)
and the permutation null as a conservative robustness check.

### Case study: an apparent topological signal of drug tolerance

We computed max-$H_1$ (top-30 PCs, $\mathbb{F}_2$) on five EGFR-mutant systems: PC9 erlotinib and
osimertinib time-series, two osimertinib PDX models, a treatment-naive patient atlas, and a
14-patient longitudinal cohort. max-$H_1$ increased with drug exposure in the osimertinib cell
line and the PDX models (osimertinib D0$\to$D14: $1.41\to3.78$; PDX YU-006 untreated$\to$residual:
$2.79\to5.46$; **Fig. 1**), the kind of monotonic "topological plasticity" trend that motivates
scTDA applications. The erlotinib series is a counter-example even at this stage — its signal is
weak and non-monotonic (D9 $1.77$, D11 $1.08$, the lowest of the series) — which we keep in view
as evidence of the statistic's instability rather than omit. Taken at face value, the cell-line
and PDX trends look like a cross-system topological signature of drug tolerance. We now subject it
to the three controls.

### Control 1 — cells do not traverse the loop

Circular coordinates from the most-persistent $H_1$ class (DSPVJ harmonic smoothing) collapse to a
single angle in every cohort: Rayleigh $R \ge 0.965$ (erlotinib D9 0.992; osimertinib D14 0.997;
PDX residual 0.965; patient-naive 0.987; patient-PD 1.000; **Fig. 2a**) — more concentrated than
the $R=0.881$ of the 12%-occupied positive control. Using the Watermelon lineage barcodes,
clonal lineages are if anything *more* angularly concentrated than size-matched random groups,
and are tighter than random at D14 — the temporal-proxy evidence that closes the snapshot caveat
(each cell is sequenced once, so traversal is inferred from clonal behaviour rather than observed
directly). The persistent loop is a real but sparse feature the population does not occupy or
traverse; max-$H_1$ does not evidence cyclic state transitions.

### Control 2 — baseline max-$H_1$ is dispersion

Against the covariance-matched Gaussian null, baseline max-$H_1$ does not exceed dispersion
(osimertinib D0/D3 ratio $\approx 0.9$; observed $>$ null in 18–31% of bootstrap subsamples;
**Fig. 1**). An excess appears only under drug (osimertinib D7/D14 observed $>$ null in 100%/98%;
PDX residual 99%). This is the control that the field-standard gene-shuffle null does not provide:
shuffling genes within cells destroys all covariance and is beaten by essentially any real data,
testing "is there any structure," not "is there structure beyond dispersion." Note the null's
failure mode is conservative for our negative conclusion: a Gaussian preserves only second
moments, so if anything it under-matches higher-moment structure and would over-call an excess —
biasing *against* the finding that the signal is dispersion.

### Control 3 — the drug excess is the cell-cycle loop

The cell cycle is intrinsically a loop and is the canonical $H_1$ confound. Localising the cells
that contribute the osimertinib D14 excess (leave-one-cluster-out on max-$H_1$), the
structure-driving population is proliferating cells: differential expression returns cell-cycle
genes (PTTG1, UBE2S, CKS1B, CENPW, CDKN3, MYBL2) and proliferation markers (MKI67, PCNA, CCNB1,
CDK1, BIRC5; 5/6), not EMT (1/19) or persister markers (0/14), and they are not low-quality cells
(median counts 10,142 vs 8,692). Decisively, in the pooled embedding with bootstrap, the drug
excess over the dispersion null does **not** survive S/G2M regression: the D14 observed/null ratio
falls from 1.18 to 0.99, and D14 $>$ D0 holds in only 31% of bootstraps after regression
(**Fig. 2b**). This exposes a methodological error in common practice: the usual cell-cycle
control reports that *raw* max-$H_1$ is preserved after regression and concludes "not cell cycle,"
but raw max-$H_1$ is dispersion-dominated; the cell-cycle contribution is visible only as a ratio
to the dispersion null.

The framework decomposes rather than uniformly debunks. In the *in-vivo* PDX system the result is
the opposite: the residual-disease excess (the strongest in the study) is **not** explained by the
cell cycle. Excluding cycling cells leaves it intact — the G1-only dispersion ratio is 2.07
(90% bootstrap CI $[1.37, 2.61]$; excess $>1$ in 100% of bootstraps), indistinguishable from a
size-matched random subset (1.95). The PDX therefore carries genuine non-Gaussian, non-cell-cycle
topological structure under drug. Crucially, this structure still fails the traversal test
($R=0.965$, Control 1): even where structure is real, the population does not occupy or traverse
it, so the "cyclic plasticity" reading remains unsupported. The identity of this in-vivo structure
(microenvironmental, mesenchymal, or otherwise) is left to future work. Across all five systems,
then: dispersion accounts for the baseline; the cell cycle accounts for the in-vitro drug excess;
the in-vivo residual structure is genuine but non-traversed — and in no system does max-$H_1$
evidence a cell-traversed cycle.

### A confound-resistant clonal statistic is under-powered on standard data

Because the above confounds act on point-cloud geometry, we asked whether a lineage-level
statistic could recover a genuine signal. We defined clonal state-memory
$M(t)=1-H_{\text{obs}}/H_{\text{null}}$, the reduction in within-clone state entropy relative to a
lineage-shuffle null preserving dispersion and cell cycle. On the Watermelon osimertinib data, a
baseline reduction in within-clone entropy is *suggested* but **not robust** (M(D0) point estimate
0.48 but 90% bootstrap CI $[-0.05, 0.29]$, crossing zero; M(D0) ranges $-0.12$ to $1.00$ across
clustering settings), and the proposed decay under drug is **not robust**: its sign flips across
clustering resolution and clone-size thresholds (positive in 3/9 settings; negative at the
most-data setting), the trajectory is non-monotonic, intermediate timepoints carry replicate and
sort confounds, and a depth-stratified null halves the apparent decay (**Fig. 3**). The root cause
is power: standard lineage-tracing data yield only 12–33 clones with $\ge 3$ cells per timepoint
(median clone size 3), too few to estimate a within-clone state distribution. This is a cautionary
result in its own right: clonal-memory statistics require clone-rich designs that current
drug-timecourse datasets rarely provide.

---

## Methods

**Datasets.** GSE134839 (PC9 erlotinib, Drop-seq), GSE150949 (PC9 osimertinib, Watermelon lineage
tracing, 10x), GSE243562 (two osimertinib PDX models, 10x), GSE131907 (treatment-naive LUAD
atlas, 10x), and the EGFR-mutant subset of Maynard et al. 2020 (Smart-seq2). Standard scanpy QC,
total-count normalisation, log1p, 3,000 HVGs, total-count and mitochondrial regression, scaling,
PCA (seed 42); groups subsampled to 1,000–1,500 cells.

**Persistent homology.** `ripser` 0.6, Vietoris–Rips to $H_1$ on top-30 PCs, $\mathbb{F}_2$;
max-$H_1=\max(\text{death}-\text{birth})$ over finite $H_1$ bars.

**Dispersion null (Control 2).** For each group, draw $n$ samples from $\mathcal{N}(\hat\mu,
\hat\Sigma)$ of the PCA cloud (15–20 draws) and recompute max-$H_1$; ratio = observed/median(null).
Bootstrap: 30 subsamples of 1,000–1,200 cells; report median, 90% CI, CV, and the fraction of
subsamples with observed $>$ null and with drug $>$ baseline. Second null: independent per-PC
permutation (preserves per-axis variance, destroys joint structure).

**Traversal test (Control 1).** `dreimac` `CircularCoords` (DSPVJ, prime 47, 400 landmarks) on the
most-persistent $H_1$ class; Rayleigh $R$ and sector coverage; validated on a noisy circle
($R=0.12$) and a 12%-occupied blob-plus-loop ($R=0.88$).

**Cell-cycle test (Control 3).** Leave-one-Leiden-cluster-out on max-$H_1$ to localise
structure-driving cells; Wilcoxon DE with gene symbols recovered from raw counts; pooled embedding
rebuilt with/without S/G2M-score regression (Tirosh sets), with the dispersion ratio bootstrapped
under both.

**Synthetic benchmark.** Four 1,200-cell scenarios (Gaussian blob; uniformly traversed loop;
sparse 270° arc on a blob; two Gaussian clusters), structure in 2-D embedded in 30-D with
isotropic noise; the full framework applied to each (`scripts/25_synthetic_benchmark.py`).

**Clonal memory.** Shared Leiden states; $M(t)=1-H_{\text{obs}}/H_{\text{null}}$ with a
within-timepoint lineage-shuffle null preserving clone sizes and the state marginal; robustness
over resolution, clone-size threshold, bootstrap, and a depth-stratified null.

**Code availability.** Reproducible from `scripts/16_*.py`–`scripts/25_*.py` (MIT license);
figures via `scripts/24_cautionary_figures.py`.

---

## Discussion

We have shown that a cross-system topological signal of drug tolerance — rising max-$H_1$ — is
largely explained by transcriptional dispersion and, in vitro, by the cell cycle; and that where
genuine non-confound structure remains (the in-vivo PDX residual), it still does not evidence
cyclic plasticity, because the population does not traverse it. The controls that establish this
are not exotic; each is the kind that, once stated, is obviously necessary. That the apparent signal nonetheless looked compelling across five
independent systems, survived published-style cell-cycle controls, and "replicated" in patient
tumours is exactly why a standard control framework matters: a confounded statistic can be
reproducible, cross-scale, and wrong.

This is a cautionary, methodological contribution rather than a biological discovery, in the
established lineage of structure-aware single-cell controls (for doublets, ambient RNA, batch, and
pseudotime). The biology of persister-state plasticity, established by lineage tracing, is not in
question; we show that static persistent-homology summaries do not capture it and instead track
nuisance variation. We also show that the obvious confound-resistant alternative — a clonal
state-memory statistic — is under-powered on standard lineage-tracing data, marking clone-rich
designs as a prerequisite for future work. The framework's validity is established
dataset-independently on synthetic ground truth; the five cancer systems are the worked case
study. We recommend that topological claims on perturbation scRNA-seq be accompanied, at minimum,
by a dispersion-matched null, a traversal test, cell-cycle regression in the dispersion-controlled
frame, and bootstrap directional testing.

**Limitations.** The cell-cycle attribution (Control 3) is demonstrated on the osimertinib cell
line; the PDX residual excess is shown to be non-cell-cycle but its positive identity is not
established here. We address max-$H_1$; the dispersion null derives from the same Vietoris–Rips
filtration and is expected to apply to other persistence summaries (persistence images, Adams et
al. 2017; landscapes, Bubenik 2015) and to higher-degree homology, though we do not test this.
Traversal is inferred from snapshot circular-coordinate concentration plus the clonal-spread
proxy, not from direct temporal observation.
