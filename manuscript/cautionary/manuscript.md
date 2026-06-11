# Structure-aware controls for topological data analysis of perturbation single-cell RNA-seq: a dispersion-matched null, a traversal test, and a cautionary case study

**Short title:** Structure-aware controls for topological analysis of scRNA-seq

**Authors:** H.-T. Lin, Y.-H. Tu

**Target journals:** Bioinformatics · GigaScience · Genome Biology (Method)

---

## Key points

- A popular topological statistic — max-$H_1$, the size of the most persistent "loop" — *appears*
  to measure cyclic drug-tolerant cell-state plasticity, and rises with drug exposure across five
  cancer systems.
- It does not survive three simple controls: **cells never traverse the loop**; the baseline signal
  is just **data spread (dispersion)**; and the drug-induced part is the **cell cycle** (which is
  itself a loop).
- We package the controls — a dispersion-matched null, a traversal test, cell-cycle regression
  judged against dispersion, and bootstrap testing — into one reusable check (**Fig. 5**), validated
  on synthetic data, and demonstrate the trap on our own analysis.
- **Takeaway:** run this check before reading any persistence-magnitude difference as biology.

---

## Abstract

As topological data analysis (TDA) enters single-cell biology, scalar summaries of persistent
homology — most simply the maximum degree-1 persistence (max-$H_1$) — are beginning to be read as
quantitative measures of "cyclic" cell-state structure or plasticity, and compared across
conditions. This is hazardous for two reasons that are individually known but not jointly
controlled in practice: Vietoris–Rips persistence magnitude grows with the spread (dispersion) and
sampling of a point cloud even in pure noise [Bobrowski2017], and the cell cycle is an intrinsic
loop in expression space [Schwabe2020, Rizvi2017]. We assemble four controls into a single,
reusable significance procedure for perturbation scRNA-seq — (i) a **covariance-matched Gaussian
null** preserving transcriptional dispersion; (ii) a **circular-coordinate traversal test** that
asks whether a detected loop is occupied before any cyclic reading; (iii) **cell-cycle regression
evaluated relative to the dispersion null** rather than on raw persistence; and (iv) **bootstrap
directional testing** — and validate it on synthetic data with known ground truth (sensitive to a
genuine traversed loop, specific against benign dispersion and non-Gaussian non-cyclic structure).
As a cautionary case study we report our own analysis of five longitudinal EGFR-mutant
lung-cancer systems, in which an apparently compelling, reproducible, cross-scale drug-associated
max-$H_1$ signal dissolves under the controls: cells never traverse the loop ($R \ge 0.965$ vs
$0.881$ for a 12%-occupied control), baseline max-$H_1$ equals the dispersion null, and the
in-vitro drug excess is the cell-cycle loop. The framework also isolates genuine structure where
it exists — a drug-emergent multiciliated program in PDX residual disease — showing it passes real
signal as well as rejecting confounded signal. We recommend these controls as standard practice as
scalar TDA readouts grow in single-cell biology, and release them as open-source code.

## Contributions and relation to prior work

The individual confounds and several remedies are known; our contribution is their **integration
into one corrected significance procedure for single-cell $H_1$**, with a validated operating
point, plus an honest worked failure.

1. A **covariance-matched Gaussian null** for max-$H_1$ on scRNA-seq, controlling for dispersion.
   The dominant published scTDA null is within-cell gene-label shuffling [Rizvi2017], which
   destroys covariance and is beaten by essentially any real data; the principled TDA-theory
   alternative is subsampling confidence sets [Fasy2014]. A dispersion-preserving (covariance-
   matched) null for persistence on single-cell expression appears not to have been used before,
   though the idea has precedent in surrogate-data nulls for time series and in weak-vs-strong-null
   critiques of geometric claims in single-cell embeddings [Kendiukhov2026].
2. A **circular-coordinate traversal test** [deSilva2011] establishing loop occupancy as a
   precondition for any "cyclic" interpretation.
3. The point — made explicit and operationalised — that **a cell-cycle control based on raw
   max-$H_1$ surviving S/G2M regression is inadequate**, because raw persistence is
   dispersion-dominated [Bobrowski2017]; the cell-cycle contribution must be judged relative to the
   dispersion null. The cell-cycle-as-loop confound itself is well established [Rizvi2017,
   Schwabe2020]; the inadequacy of the raw-persistence control, to our knowledge, is not.
4. A **synthetic benchmark** characterising sensitivity and specificity, and a five-system
   **cautionary case study** — our own — in which a confounded signal is reproducible, cross-scale,
   and wrong, demonstrating that the trap is real for careful practitioners.

---

## Introduction

Topological data analysis is entering single-cell biology, with reviews now cataloguing
applications of persistent homology and Mapper to scRNA-seq [HernandezLemus2025]. Its appeal is
the ability to detect closed, non-tree-like structure (loops; degree-1 homology, $H_1$) that
trajectory-inference tools, designed for one-directional progressions, do not quantify. As the
field matures, a natural step is to summarise topology with a scalar — max-$H_1$, total
persistence, persistence entropy — and to compare that scalar across conditions, reading an
increase as more cyclic structure or plasticity. In cancer drug tolerance this is tempting:
lineage-tracing shows persister cells reversibly transition among states, a behaviour readily
called "cyclic."

We argue this step needs controls that current practice does not jointly apply, and we provide
them. Two confounds are individually documented. First, Vietoris–Rips persistence is sensitive to
the spread and sampling density of the point cloud: the maximally persistent cycle of pure random
points grows with sample size [Bobrowski2017], raw persistence magnitude is metric- and
scale-dependent (motivating stabilised summaries such as persistence images and landscapes
[Adams2017, Bubenik2015]) and outlier-sensitive (motivating density-robust filtrations
[Anai2020]). An increase in transcriptional dispersion under perturbation can therefore raise
max-$H_1$ with no change in genuine topology. Second, the cell cycle traces a circle in expression
space [Schwabe2020] and is the canonical generator of an $H_1$ loop in scRNA-seq; the
field-defining scTDA method states explicitly that "the cell cycle will give rise to periodic
structures in the expression space" [Rizvi2017]. The nulls in common use do not control for the
first (gene-label shuffling [Rizvi2017] destroys covariance, including dispersion), and a recent
review notes that single-cell TDA "may not provide rigorous statistical assessments" and is
"sensitive to preprocessing" without resolving it [HernandezLemus2025].

This is a preventive contribution, not a claim that a widespread error is being committed: scalar
persistence-vs-condition comparison on scRNA-seq is still uncommon (the established scTDA methods
use Mapper graphs and per-cell features [Rizvi2017, Nicolau2011]), and careful applications that
control sample size and use diagram distances exist [Mukherjee2022]. Our point is that as the
scalar-readout step becomes attractive, the controls below should accompany it — and that the trap
is real, because we fell into it ourselves. We present our own five-system analysis as the
cautionary example.

---

## Results

### The control framework

The framework takes a cell-by-PC embedding per condition and reports three quantities plus a
decision rule (**Fig. 5**). In plain terms: persistence magnitude measures how big the most
robust "loop" in the data is, but a loop can be faked by cells simply being more spread out
(*dispersion*) or by the cell cycle (which traces a circle in expression space), so we ask three
questions in turn — is there more loop than spread alone would give? do cells actually go around
it? and is what's left just the cell cycle? (1) *Dispersion test:* the ratio of observed max-$H_1$ to that of a multivariate
Gaussian matched to the embedding's mean and covariance; a ratio near 1 means the statistic is
explained by dispersion, a ratio robustly $>1$ indicates structure beyond dispersion.
(2) *Traversal test:* the Rayleigh concentration $R$ of circular coordinates [deSilva2011] from the
most-persistent $H_1$ class; $R$ near 1 means the population sits at one angle (loop not occupied).
(3) *Cell-cycle test:* the dispersion ratio recomputed after S/G2M-score regression
[Tirosh2016]. All three are wrapped in bootstrap subsampling [Fasy2014]; claims rest on the
fraction of subsamples in which an ordering holds. Decision rule: *structure beyond dispersion* if
the Gaussian ratio is robustly $>1$; *a traversed loop* only if additionally $R<0.6$;
*cell-cycle driven* if the excess does not survive S/G2M regression.

### The framework is sensitive and specific on synthetic ground truth

On point clouds with known answers (n=1,200 each, structure in two dimensions embedded in 30-D with
isotropic noise; **Fig. 4**): a Gaussian blob gives ratio 0.93, $R=0.99$ (correctly *no
structure*); a uniformly traversed loop gives ratio 8.35, $R=0.42$ (correctly a *traversed loop*;
sensitivity); a sparse non-closed arc gives ratio 1.54, $R=0.96$ (*structure, not traversed*); two
separated Gaussian clusters — non-Gaussian but non-cyclic — give ratio 0.93 (*no structure*;
specificity). The Gaussian null does not false-positive on benign non-Gaussianity, the key
objection to a covariance-matched null. A second non-parametric dispersion-preserving null
(per-PC permutation) agrees and is more conservative.

### Case study: an apparent topological signal of drug tolerance

We computed max-$H_1$ (top-30 PCs, $\mathbb{F}_2$) on five EGFR-mutant systems. max-$H_1$ increased
with drug in the osimertinib cell line and PDX models (osimertinib D0$\to$D14 $1.41\to3.78$; PDX
YU-006 untreated$\to$residual $2.79\to5.46$; **Fig. 1**). The erlotinib series is already a
counter-example — weak and non-monotonic (D9 1.77, D11 1.08, the lowest of the series).

### Control 1 — cells do not traverse the loop

Circular coordinates collapse to one angle in every cohort: Rayleigh $R \ge 0.965$ (**Fig. 2a**) —
more concentrated than the $R=0.881$ of a 12%-occupied positive control. Clonal lineages (Watermelon
barcodes) are if anything *more* concentrated than random. The loop is a sparse feature the
population does not occupy or traverse.

### Control 2 — baseline max-$H_1$ is dispersion

Against the covariance-matched Gaussian null, baseline max-$H_1$ does not exceed dispersion
(osimertinib D0/D3 ratio $\approx 0.9$; observed $>$ null in 18–31% of subsamples; **Fig. 1**). An
excess appears only under drug (osimertinib D7/D14 in 100%/98%; PDX residual 99%). Because a
Gaussian preserves only second moments, its failure mode is conservative for a negative
conclusion.

### Control 3 — the in-vitro excess is the cell cycle; the in-vivo excess is not

The osimertinib D14 excess is contributed by proliferating cells (DE: PTTG1, UBE2S, CKS1B, MYBL2;
proliferation markers MKI67/PCNA/CCNB1/CDK1/BIRC5, 5/6) and does not survive S/G2M regression: the
D14 ratio falls 1.18$\to$0.99 and D14 $>$ D0 holds in only 31% of bootstraps (**Fig. 2b**). The
common control — *raw* max-$H_1$ preserved after regression — would have wrongly concluded "not cell
cycle," because raw max-$H_1$ is dispersion-dominated [Bobrowski2017]. The framework decomposes
rather than uniformly debunks: in the in-vivo PDX the residual excess is *not* cell cycle (G1-only
ratio 2.07, 90% CI $[1.37,2.61]$, $>1$ in 100% of bootstraps), yet still fails the traversal test
($R=0.965$). Across all five systems: dispersion accounts for the baseline; the cell cycle for the
in-vitro excess; the in-vivo residual is genuine but non-traversed — and in no system does
max-$H_1$ evidence a cell-traversed cycle.

**Summary of the case study.** What each control returns per system:

| System | Dispersion test | Traversal ($R$) | Cell-cycle test | Verdict |
|--------|-----------------|-----------------|-----------------|---------|
| Osimertinib cell line | baseline $\approx$1, drug 1.4 | 0.997 (not traversed) | excess **is** cell cycle | confounded |
| PDX (YU-006) | drug 1.95 | 0.965 (not traversed) | excess **not** cell cycle | genuine but non-traversed |
| Erlotinib cell line | $\approx$1 (no excess) | — | — | no signal |
| Patient tumours (Maynard/Kim) | excess at baseline too | high | n/a | heterogeneity, not drug |

### A confound-resistant clonal statistic is under-powered

A clonal state-memory statistic $M(t)=1-H_{\text{obs}}/H_{\text{null}}$ against a lineage-shuffle
null is non-robust on the Watermelon data: a baseline reduction is suggested but not robust (M(D0)
point estimate 0.48, 90% CI $[-0.05,0.29]$), and its decay under drug flips sign across clustering
settings (positive in 3/9), is non-monotonic, and halves under a depth-stratified null
(**Fig. 3**). The cause is power: only 12–33 clones with $\ge 3$ cells per timepoint.

### Beyond rejection: the framework isolates a genuine, interpretable program

A control framework earns trust only if it passes real structure (cf. the synthetic true positive).
The PDX residual is the worked example: its excess, the one signal surviving every control, is
driven by a coherent multiciliated airway-epithelial program (PIFO, RSPH1, CFAP45/126/157,
CAPS/CAPSL, CETN2, TPPP3; secretory AGR2/AGR3), drug-emergent in YU-006 (0.1% of untreated cells
$\to$ 9.6% under residual disease; $\approx$6 $\to$ $\approx$187 cells in absolute terms while the
tumour contracts threefold), consistent with the differentiated lineage states EGFR-TKI persisters
adopt. Two bounds: the program is model-specific (absent in YU-003), and copy-number inference
could not resolve tumour-derived (transdifferentiation) vs an expanded non-malignant population
without a matched-normal reference. This real structure still fails the traversal test.

---

## Methods

**Datasets.** GSE134839 (PC9 erlotinib, Drop-seq), GSE150949 (PC9 osimertinib, Watermelon, 10x),
GSE243562 (PDX, 10x), GSE131907 (LUAD atlas, 10x), Maynard et al. 2020 EGFR subset (Smart-seq2);
standard scanpy QC/normalisation/HVG/PCA (seed 42), subsampled to 1,000–1,500 cells.
**Persistent homology.** `ripser` 0.6, Vietoris–Rips to $H_1$ on top-30 PCs, $\mathbb{F}_2$.
**Dispersion null.** draws from $\mathcal{N}(\hat\mu,\hat\Sigma)$ of the PCA cloud; ratio $=$
observed/median(null); 30 bootstrap subsamples, 90% CIs; second null by per-PC permutation.
**Traversal test.** `dreimac` `CircularCoords` (DSPVJ [deSilva2011]), Rayleigh $R$; validated on a
noisy circle ($R=0.12$) and a 12%-occupied loop ($R=0.88$). **Cell-cycle test.**
leave-one-cluster-out, Wilcoxon DE, S/G2M regression [Tirosh2016] with the dispersion ratio
bootstrapped. **Synthetic benchmark.** four 1,200-cell scenarios. **In-vivo identity.** DE on raw
counts; ciliated signature scored across PDX models/conditions; CNV by `infercnvpy`.
**Code.** `scripts/16_*.py`–`scripts/27_*.py` (MIT).

---

## Discussion

A cross-system topological signal of drug tolerance — rising max-$H_1$ — is largely explained by
transcriptional dispersion and, in vitro, by the cell cycle; where genuine non-confound structure
remains, it still does not evidence cyclic plasticity, because the population does not traverse it.
That the apparent signal looked compelling across five systems, survived a published-style
cell-cycle control, and "replicated" in patient tumours is exactly why these controls matter: a
confounded statistic can be reproducible, cross-scale, and wrong.

Our contribution is preventive and integrative. Each confound has prior literature — density- and
sample-size-dependence of persistence [Bobrowski2017], cell-cycle loops [Schwabe2020, Rizvi2017],
preprocessing sensitivity and the absence of significance testing in single-cell TDA
[HernandezLemus2025], and, closest to our thesis, the demonstration in single-cell foundation-model
embeddings that geometric structure significant under weak (shuffle) nulls vanishes under stronger
nulls [Kendiukhov2026]. We differ from that work in target and null (raw-expression $H_1$ with a
dispersion-matched Gaussian null, vs embedding geometry with rewiring nulls) and contribute the
*integration* of dispersion null, traversal test, and dispersion-relative cell-cycle control into
one procedure for single-cell persistence, validated on ground truth. The framework also isolates
genuine structure, surfacing a drug-emergent multiciliated program in residual disease; that lead,
and the question of whether it is tumour transdifferentiation conferring tolerance, is a
prospective, experiment-bearing programme rather than a reanalysis.

**Outlook.** Establishing the in-vivo program as a discovery would require resolving cell-of-origin
(matched-normal CNV or EGFR-driver genotyping), reproduction in further models or patient
residual-disease specimens, and functional testing of whether the differentiated state confers
tolerance.

**Limitations.** The cell-cycle attribution is demonstrated on the osimertinib cell line; the PDX
residual excess is shown non-cell-cycle but its cell-of-origin is unresolved. We address max-$H_1$;
the dispersion null derives from the same filtration and should apply to other persistence
summaries [Adams2017, Bubenik2015] and higher-degree homology, untested here. Traversal is inferred
from snapshot circular-coordinate concentration plus the clonal-spread proxy.

---

## References

[Rizvi2017] Rizvi AH, Cámara PG, Kandror EK, et al. Single-cell topological RNA-seq analysis reveals insights into cellular differentiation and development. *Nat Biotechnol* 2017;35:551–560. doi:10.1038/nbt.3854.

[Nicolau2011] Nicolau M, Levine AJ, Carlsson G. Topology based data analysis identifies a subgroup of breast cancers with a unique mutational profile and excellent survival. *PNAS* 2011;108:7265–7270.

[Bobrowski2017] Bobrowski O, Kahle M, Skraba P. Maximally persistent cycles in random geometric complexes. *Ann Appl Probab* 2017;27:2032–2060.

[Fasy2014] Fasy BT, Lecci F, Rinaldo A, Wasserman L, Balakrishnan S, Singh A. Confidence sets for persistence diagrams. *Ann Statist* 2014;42:2301–2339.

[Schwabe2020] Schwabe D, Formichetti S, Junker JP, Falcke M, Rajewsky N. The transcriptome dynamics of single cells during the cell cycle. *Mol Syst Biol* 2020;16:e9946.

[Adams2017] Adams H, Emerson T, Kirby M, et al. Persistence images: a stable vector representation of persistent homology. *JMLR* 2017;18:1–35.

[Bubenik2015] Bubenik P. Statistical topological data analysis using persistence landscapes. *JMLR* 2015;16:77–102.

[Anai2020] Anai H, Chazal F, Glisse M, et al. DTM-based filtrations. In *Topological Data Analysis* (Abel Symposia 15), Springer, 2020 (cf. Chazal et al., Robust Topological Inference, *JMLR* 2017;18:1–40).

[HernandezLemus2025] Hernández-Lemus E. Topological data analysis in single cell biology. *Front Immunol* 2025;16:1615278.

[Kendiukhov2026] Kendiukhov I. What topological and geometric structure do biological foundation models learn? Evidence from 141 hypotheses. *arXiv* 2026;2602.22289.

[deSilva2011] de Silva V, Morozov D, Vejdemo-Johansson M. Persistent cohomology and circular coordinates. *Discrete Comput Geom* 2011;45:737–759.

[Tirosh2016] Tirosh I, Izar B, Prakadan SM, et al. Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq. *Science* 2016;352:189–196.

[Mukherjee2022] Mukherjee S, Wethington D, Dey TK, Das J. Determining clinically relevant features in cytometry data using persistent homology. *PLoS Comput Biol* 2022;18:e1009931.
