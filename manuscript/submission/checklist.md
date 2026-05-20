# Submission checklist — *Nature Communications*

Tracking the items the editorial team will check on receipt. Verify each
before clicking "Submit" at <https://submission.springernature.com/>.

## Manuscript file

- [ ] **Title** ≤ 15 words, no abbreviations
  - current: "Drug tolerance takes shape: persistent homology reveals cyclic cell-state plasticity in EGFR-mutant lung cancer" — **14 words ✓**
- [x] **Running title** ≤ 50 characters — *"Topology of drug-tolerant lung cancer plasticity"* (47 chars) — enter in portal metadata; nature.cls does not render it
- [ ] **Abstract** ~ 150 words, unstructured, no citations — *currently structured (Background/Results/Conclusions); de-structure for NC*
- [ ] **Main text + Methods ≤ 5,000 words**
  - current `texcount`: **7,835 words** — `WORK REMAINING`: ~2,800 words to trim
- [ ] **Display items ≤ 10 figures + tables** combined
  - current: 9 figures + 5 main tables = **14 items** — `WORK REMAINING`: consolidate or move tables to SI
- [ ] **References ≤ 60** numbered (Nature style)
  - current bibliography: see `references.bib` — verify count and trim if needed
- [ ] **Figure resolution** ≥ 300 dpi photos / ≥ 1200 dpi line art / ≥ 600 dpi mixed
  - all figures produced at 300 dpi ✓ (matplotlib `savefig.dpi=300`); confirm vector PDFs preferred for line art
- [ ] **Figure file format:** PDF (vector) or TIFF/EPS for raster; per-figure files
- [ ] **Figure 1 AI-image compliance:** vector matplotlib, no AI-generated submission asset ✓

## Author / affiliation

- [ ] First author: **Hsieh-Ting Lin, M.D.** (KFSYSCC) ✓
- [ ] Corresponding author: **Yueh-Hua Tu, Ph.D.** (CCSB NTU) ✓
- [x] Corresponding email entered in author block — `yuehhuatu@ntu.edu.tw`
- [ ] ORCID iDs for both authors entered in submission form
- [ ] Author contribution statement (CRediT) attached — see `credit_taxonomy.md` ✓

## Declarations

- [ ] Competing interests — none ✓
- [ ] Funding — **`<to be filled in by Y-H Tu>`**
- [ ] Ethics — public datasets only; no new human/animal work ✓
- [ ] Data availability — GEO accessions listed; Zenodo deposit to follow ✓
- [ ] Code availability — GitHub MIT-licensed ✓; Zenodo software DOI to follow
- [ ] AI-tool usage statement included ✓

## Cover letter

- [ ] Cover letter attached — see `cover_letter.md` ✓
- [ ] Suggested reviewers (5–10) — 9 listed in cover letter ✓
- [ ] Excluded reviewers — none
- [ ] Not under consideration elsewhere — confirmed ✓
- [ ] Preprint disclosed — bioRxiv (DOI to be added at upload)

## Pre-submission integrity

- [ ] Spell-check / grammar pass — pending Vale lint run
- [ ] All `\ref{}` resolve, no `??` in PDF — verified ✓
- [ ] All figures render with current `\includegraphics{...}` paths — verified ✓
- [ ] `references.bib` lints cleanly — verify before submit
- [ ] Inference-audit critiques addressed — see `manuscript/audit/inference_audit.argdown`:
  - [ ] Replace "isometry-invariant" with metric-aware wording
  - [ ] Cite Oren 2021 for the dynamic-cycling claim in C6 / abstract
  - [ ] Verify the "6/6 p<0.05, 5/6 p<0.01" claim against `tab:validation`
  - [ ] Soften EMT conclusion to "include EMT programs" or note filter dependence
  - [ ] Mark mechanistic speculation explicitly with "We speculate that..."
- [ ] Word count trimmed to ≤ 5,000 (D-prose-polish workstream)

## Zenodo / GitHub release pre-submission

- [ ] Tag `v1.0.0` release on GitHub
- [ ] Zenodo–GitHub integration triggers software DOI
- [ ] Update `data availability` and `code availability` blocks with both DOIs
- [ ] Final manuscript PDF rebuild (`make pdf`) after DOI insertion

## Format compliance

- [ ] LaTeX source migrated to Springer Nature `sn-jnl.cls` (option `nature`)
  - working: `manuscript/preprint.tex` (article class, bioRxiv build)
  - target: `manuscript/preprint_ncomms.tex` (sn-jnl class) — see `template_setup.md`
- [x] Line numbering ON for review version (Nature convention) — `\usepackage{lineno}` + `\linenumbers` added to preprint_ncomms.tex preamble; remove for camera-ready
- [ ] Reference style: numbered, Nature style (`unsrtnat` or sn-nature.bst)
