---
title: "Phase 0 result — feasibility scan for the topology audit"
date: 2026-06-12
decision: "NO-GO on the audit thesis; N_auditable far below gate"
gate: "GO only if N_auditable >= 6-8 vulnerable single-cell findings with public data"
tags: [phase0, go-no-go, audit, decision]
---

# Phase 0 result: NO-GO on the audit

Four independent literature+data scans (PH-magnitude, Mapper, cyclic-trajectory, recent/review
sweep). They converge.

## Count

**Genuinely vulnerable + auditable single-cell "topological loop as biology" findings: ~1–3
(1 clean, 2 weak/borderline) — far below the gate of 6–8.**

| Candidate | Claim | Vulnerable? | Auditable? | Note |
|-----------|-------|-------------|-----------|------|
| **TopGen** (Flores-Bautista & Thomson 2023, bioRxiv) | H1 loops = stem/transdiff/"circuits" in >10 dev atlases | yes (no dispersion/CC control) | **yes** | strongest target — but a **preprint**, developmental not cancer |
| **TopoCytometry** (Mukherjee 2022, PLoS CB) | H0/H1 disease-state topology, COVID vs healthy | partial | maybe | cytometry not scRNA; actually controls sample size (careful) |
| **Wang 2019** (PSB, Mapper melanoma) | Mapper "continuity" in melanoma | weak | maybe | claim is continuity, not a discrete loop |
| Rizvi 2017 (scTDA, GSE94883) | β1 loops in differentiation | **NO — positive control** | n/a | **DID** the cell-cycle representation + PH null; the loop IS the cell cycle. The exemplar, not a target. |
| Oren 2021, Aissa 2021, Chauvistré 2022 | "cycling/reversible" persisters | excluded | — | "cycling" = cell-cycle re-entry, or did the CC control (Aissa). Not loop-in-state-space. |
| Nicolau 2011, Rabadan 2020 | Mapper flares/structure | — | — | **bulk**, not single-cell. Out of scope. |

## Why NO-GO (the convergent findings)

1. **The vulnerable practice barely exists in published single-cell work.** Genuine "uncontrolled
   topological-loop-as-biology" single-cell claims are rare (~1 clean, a preprint).
2. **The field-defining method did the controls right.** Rizvi 2017 (scTDA) built a
   cell-cycle-gene-only representation, found the loop *was* the cell cycle, and ran a PH null —
   it is the positive control, not a victim. An audit would mostly vindicate the field.
3. **The "cycling persister" keyword is a trap** — it means cell-cycle re-entry/proliferation
   (Oren), or the paper already controlled cell cycle (Aissa). Not closed loops in state space.
4. So "many published topological findings are artifacts" **cannot be supported** — there aren't
   many to audit, and auditing papers that controlled correctly (or never claimed a loop) would
   itself be a straw man.

## Decision

**NO-GO on the audit thesis.** Do not pivot. The gate did its job: ~1–2 weeks of scanning saved
~3–6 months on a project whose target set does not exist.

## Two genuinely useful by-products (fold into the cautionary paper)

1. **Empirical backing for the preventive framing.** We can now state, with evidence from a
   systematic scan, that uncontrolled topological-loop claims in single-cell are rare and the
   canonical method (Rizvi 2017) controls correctly — so our contribution is *preventive*, not a
   correction of widespread malpractice. This hardens the de-strawmanning with data.
2. **Our niche is confirmed unoccupied.** All four scans independently found no prior single-cell
   TDA paper on cyclic plasticity in EGFR/PC9 drug tolerance — the cautionary paper is first in
   its niche and not pre-empted.

## Recommendation

Ship the cautionary/methods paper (preventive framing, already de-strawmanned), and add one
sentence + a supplementary note citing this scan as the empirical basis for "preventive, not
corrective." Drop the audit.
