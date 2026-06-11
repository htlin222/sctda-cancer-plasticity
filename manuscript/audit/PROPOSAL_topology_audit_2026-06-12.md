---
title: "Proposal — A reproducibility audit of topological 'structure' claims in single-cell RNA-seq"
date: 2026-06-12
status: proposal (pre go/no-go); pilot already complete
relation: extends the cautionary/methods paper (manuscript/cautionary/) from preventive to audit
author: H.-T. Lin
tags: [proposal, audit, reproducibility, tda, go-no-go]
---

# Proposal: Do published topological "structures" in single-cell RNA-seq survive structure-aware controls?

**Working title:** *Reproducible but not real? A structure-aware audit of topological cell-state
findings in single-cell RNA-seq.*

## 1. One-paragraph thesis

Single-cell studies increasingly report topological cell-state "structure" — loops, cycles,
non-tree manifolds — detected by persistent homology (persistence magnitude / diagrams), Mapper, or
cyclic-trajectory methods. We have shown, in a five-system pilot, that such a signal can be
reproducible across cell line, PDX, and patient tumours and *still* be an artefact of two confounds
(transcriptional dispersion and the cell cycle). This proposal asks the obvious next question
**empirically and at scale**: of the published single-cell topological-structure claims that used a
vulnerable method and have public data, **what fraction survive a structure-aware control
framework?** The answer — whether 10% or 90% — is itself the contribution.

## 2. Why this is interesting (and why now)

- It elevates the work from *preventive* ("controls should be standard") to *diagnostic* ("here is
  how often the field's topological claims hold up"). The reproducibility-audit genre (cf. doublet,
  ambient-RNA, batch-effect critiques) has had outsized impact **when the audited practice is real
  and the finding is quantitative**.
- We already have the instrument: a validated control framework (dispersion-matched Gaussian null,
  circular-coordinate traversal test, dispersion-relative cell-cycle control, bootstrap), with a
  synthetic sensitivity/specificity benchmark, and a worked pilot (our own analysis).
- It is honest by construction: we audit *our own* finding first, as the index case.

## 3. The pilot (already done) = proof of concept

Our EGFR-mutant five-system analysis: max-$H_1$ rose monotonically with drug, survived a
published-style cell-cycle control, and "replicated" in patient tumours — yet under the framework,
cells never traverse the loop ($R\ge0.965$), baseline equals the dispersion null, and the
drug-induced excess is the cell cycle. One genuine non-confound signal survived (a drug-emergent
ciliated program in PDX) but is itself non-traversed. This is the index case the audit generalises.

## 4. Central questions

1. **Prevalence:** how many published single-cell topological-structure claims used a vulnerable
   method (persistence magnitude across conditions; Mapper/PH loops interpreted as biology) without
   dispersion and cell-cycle controls?
2. **Survival:** of those with re-analysable public data, what fraction survive the framework, and
   into which class do failures fall (dispersion / cell cycle / non-traversed)?
3. **Predictors:** what distinguishes surviving claims (occupied loops, lineage support, validation)
   from failing ones — i.e. what should reviewers ask for?

## 5. Scope: inclusion / exclusion

**Include:** single-cell (or comparable high-dimensional omics) studies that present a
topological/cyclic/loop/non-tree structure **as a biological finding**, via persistent homology
(magnitude or diagram comparison), Mapper loop/flare claims, or persistence-based cyclic trajectory.

**Exclude:** correct cell-cycle *detection* (the loop is the cell cycle, by design); pure
clustering / dimensionality-reduction / classification uses with no structural biological claim;
studies with no recoverable data.

## 6. Approach (phased, gated)

**Phase 0 — Feasibility scan (the go/no-go gate; ~1–2 weeks).**
Systematic literature enumeration (PubMed / bioRxiv / citation graphs of Rizvi 2017, Nicolau 2011,
giotto-tda / kepler-mapper / ripser usage) → for each candidate, record: the structural claim, the
method, whether controls were applied, and **public data availability + re-analysability**. Output:
a count `N_auditable` and a triaged target list.
→ **GO if `N_auditable` ≥ 6–8 with usable public data; otherwise NO-GO** (keep the preventive paper).

**Phase 1 — Audit (per finding; ~3–5 months).**
For each auditable claim: obtain data, reproduce the reported structure (confirm we can recover what
they saw), then run the framework and classify the result. Pre-registered protocol and fixed
thresholds; blinded where possible.

**Phase 2 — Synthesis.**
Prevalence estimate with uncertainty; survival classification; predictors of survival; a reviewer
checklist. Release the per-finding reproducibility reports and the toolkit.

## 7. Deliverables

1. The audit paper (prevalence + survival + recommendations).
2. The open-source control framework (already built; `scripts/16–28`).
3. A public, per-finding reproducibility appendix (data, code, verdict) — the credibility backbone.

## 8. Risks and mitigations (honest)

| Risk | Severity | Mitigation |
|---|---|---|
| **Target set too thin** (our prior scan: the narrow practice is rare) | High | Phase 0 gate decides before any commitment; widen scope to Mapper/cyclic-trajectory claims, not just max-$H_1$ |
| **Confrontational** (naming published work as artefact) | High | Frame as "cannot be distinguished from confound X," not "is wrong"; pre-register; offer authors right-of-reply; audit our own work first |
| **Data unavailable / un-reproducible** | Medium | Triage in Phase 0; report what could not be audited transparently (no silent dropping) |
| **Novelty cap** (Kendiukhov 2026 etc. adjacent) | Medium | Position as the single-cell-specific, at-scale empirical audit, not the conceptual first |
| **Reproducing others' pipelines is slow** | Medium | Bound to N findings; standardise the framework application |

## 9. Go/no-go gates

- **Gate A (after Phase 0):** proceed only if `N_auditable` ≥ 6–8 with public data. Else publish the
  preventive cautionary paper as-is.
- **Gate B (after first ~3 audits):** if a non-trivial fraction fail, continue; if all survive
  cleanly, pivot the paper to "topological claims are mostly robust — here is why" (still publishable,
  opposite headline).

## 10. Effort, venue, expected value

- **Effort:** Phase 0 ~1–2 weeks (low cost); full audit ~3–6 months.
- **Venue ceiling:** if the finding is striking and bulletproof — Genome Biology, eLife, Nature
  Methods (Analysis/Brief Comm). If modest — Bioinformatics / GigaScience.
- **EV:** higher ceiling, higher variance than the current paper. The cautionary/methods paper is the
  **floor and should be submitted regardless**; the audit is an **upside bet, fully gated on Phase 0**.

## 11. Relationship to the current paper

The current cautionary/methods paper (`manuscript/cautionary/`) is self-contained and submission-
ready. This proposal does **not** replace it; it either (a) absorbs it as the index case of a larger
audit (if Phase 0 passes), or (b) leaves it untouched (if Phase 0 fails). Recommended: **submit the
cautionary paper now; run Phase 0 in parallel; decide the audit on the count.**

## 12. Immediate next step

Run **Phase 0** — the feasibility scan producing `N_auditable`. This is the single number that turns
"is the audit worth it?" from opinion into a decision. ~1–2 weeks of literature + data triage; can be
started immediately and does not block the cautionary-paper submission.
