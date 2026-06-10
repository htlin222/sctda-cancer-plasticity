# Cover letter

**To:** The Editors, *Bioinformatics* (Discovery Notes / Applications)
**Re:** "A control framework for topological data analysis of perturbation single-cell RNA-seq: dispersion-matched nulls, traversal tests, and a drug-tolerance case study"

Dear Editors,

We submit a methods contribution addressing a confound that is becoming consequential as
topological data analysis (TDA) spreads through single-cell genomics. The maximum degree-1
persistence (max-$H_1$) is increasingly used as a scalar readout of "cyclic" cell-state
plasticity, but it conflates genuine topology with two mundane confounds — transcriptional
dispersion and the cell cycle — and the nulls in common use (within-cell gene shuffles) do not
control for either.

We provide a reusable four-part control framework: a covariance-matched Gaussian null that
separates real structure from dispersion; a circular-coordinate traversal test that asks whether a
detected loop is actually occupied before any cyclic claim; cell-cycle regression evaluated within
the dispersion-controlled frame (we show the common control — raw max-$H_1$ preserved after
regression — is inadequate because raw persistence is dispersion-dominated); and bootstrap
directional testing. On synthetic data with known ground truth the framework is sensitive (it
passes a genuine traversed loop) and specific (it rejects benign dispersion and non-Gaussian
non-cyclic structure, the key objection to a covariance-matched null).

As a worked case study we apply it to five longitudinal EGFR-mutant lung-cancer systems, where an
apparently compelling, cross-system, drug-associated max-$H_1$ signal — reproducible and present
from cell line to PDX to patient tumour — fails every control: cells do not traverse the loop, the
baseline signal equals a dispersion null, and the drug-induced excess is contributed by
proliferating cells and does not survive cell-cycle regression. The lesson is general: a
confounded topological statistic can be reproducible, cross-scale, and wrong.

We believe this fits *Bioinformatics* as a negative result paired with a reusable correction,
useful to the growing community applying TDA to single-cell data. All analyses are reproducible
from open-source code (MIT license), the framework runs unmodified on Drop-seq, 10x, and
Smart-seq2 data, and we declare no competing interests. The work has not been published elsewhere
and is not under consideration by another journal.

Sincerely,
H.-T. Lin and Y.-H. Tu
