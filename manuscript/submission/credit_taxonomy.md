# Author contributions (CRediT taxonomy)

CRediT roles per <https://credit.niso.org/>. Verbs listed in the order of involvement.

## Hsieh-Ting Lin (KFSYSCC, first author)

- **Conceptualisation** — formulation of the topological-plasticity hypothesis for EGFR-TKI tolerance.
- **Methodology** — design of the persistent-homology / Mapper pipeline; choice of statistical controls (permutation null, Tirosh and stringent cell-cycle ablations, Wasserstein comparison, gene-attribution scoring).
- **Software** — implementation of the `sctda-cancer-plasticity` pipeline; all 15 numbered modules; unit-test suite.
- **Formal analysis** — execution of persistent-homology, Mapper, permutation, Harmony, and enrichment analyses on all five datasets.
- **Investigation** — dataset curation (GSE134839, GSE150949, GSE243562, GSE131907, Maynard 2020); QC; sensitivity sweeps.
- **Data curation** — preprocessing of all AnnData objects; persistence diagrams; Mapper graphs.
- **Writing — original draft** — Background, Results, Methods, Discussion, Conclusions, all figure legends.
- **Visualisation** — all main and supplementary figures (Figs 1–9 and SI figures).

## Yueh-Hua Tu (CCSB NTU, corresponding author)

- **Conceptualisation** — research-question framing; positioning within the cancer-plasticity literature.
- **Methodology** — scientific design review; validation strategy across PC9 / PDX / patient cohorts.
- **Supervision** — analysis design, scientific direction, manuscript scope.
- **Writing — review & editing** — substantive revisions across all sections; framing for *Nature Communications* audience.
- **Project administration** — submission management; reviewer-correspondence coordination.
- **Funding acquisition** — `<to be filled in by Y-H Tu>`.

## Notes

- Both authors have read and approved the final manuscript.
- No professional or AI-generated text was used in the main manuscript prose without author review and editing. The conceptual Figure 1 was sketched via an AI tool (OpenAI `gpt-image-1`) for compositional reference only and was then re-implemented from scratch as a vector matplotlib figure — the submitted Figure 1 is wholly author-authored vector art.
- Software development assistance (Claude Code by Anthropic) was used for code-pipeline scaffolding, refactoring, and figure-script implementation; all scientific decisions, statistical analyses, and prose are the authors' own.
