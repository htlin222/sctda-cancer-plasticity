# Submission checklist — Bioinformatics (Oxford)

**Manuscript:** A control framework for topological data analysis of perturbation single-cell
RNA-seq. **Files:** `manuscript.md` (content), `preprint_cautionary.pdf` (compiled 6-page),
`preprint_cautionary.tex` (source), `cover_letter.md`, figures in `../../figures/cautionary/`.

## Venue facts (verified 2026-06)
- Bioinformatics is **fully open access** (no longer hybrid) — an APC is mandatory.
- APC ≈ **US$3,600–3,800** via OUP's live SciPris calculator (research-article rate; the public
  example figure is $3,625). **−15%** for ISCB members; LMIC + hardship waivers exist.
- **Check NTU's OUP Read & Publish agreement first** — if covered, APC = $0.
- Q1 journal (Mathematical & Computational Biology; Biochemical Research Methods).
- Suggested article type: **Original Paper** (or **Discovery Note** if trimmed) under the
  *Bioinformatics of Disease* / *Gene expression* subject area.

## Pre-submission to-do (author)
- [ ] Confirm OUP Read & Publish coverage via NTU library (decides APC cost).
- [ ] Fill author affiliations, ORCIDs, corresponding-author email, funding statement.
- [ ] Convert `manuscript.md` to the journal's required format if needed (the `.tex` preprint is
      submission-grade; Bioinformatics accepts PDF + LaTeX).
- [ ] Confirm data availability statement (all GEO accessions are public; Maynard via authors).
- [ ] Confirm code availability (GitHub repo public; tag a release for the DOI).
- [ ] Competing interests: none. Author contributions (CRediT).
- [ ] Suggest 3–4 reviewers (TDA-for-single-cell + scRNA methods).

## Strengths to foreground in the cover letter (already drafted)
1. Reusable control framework (dispersion null + traversal test + cc-in-dispersion-frame + bootstrap).
2. The methodological correction: raw-max-H1-after-regression is an inadequate cell-cycle control.
3. Synthetic ground-truth validation (sensitive + specific).
4. Worked case study: a reproducible, cross-scale, *wrong* signal — and the framework also
   isolates a genuine, interpretable in-vivo program (proves it passes real structure too).

## Known reviewer questions & where they are answered
- "Is the Gaussian null too strong?" → synthetic two-cluster specificity test + second per-PC null.
- "Single-system cell-cycle claim?" → bounded to the cell line; PDX shown non-cell-cycle.
- "Just a null reanalysis?" → toolkit-forward framing + synthetic benchmark + the positive
  in-vivo isolation result.
- "Generality beyond max-H1 / EGFR?" → bounded in Limitations; synthetic validity is method-level.

## Honest expectation
Q1, realistic acceptance after normal revision. Not a high-IF/headline venue — that path needs
the experiment-bearing follow-up flagged in the Outlook, not reanalysis.
