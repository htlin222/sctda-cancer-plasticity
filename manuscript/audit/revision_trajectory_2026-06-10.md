---
title: "Revision Trajectory Assessment — Converging or Oscillating?"
date: 2026-06-10
scope: 58 commits, 2026-04-14 → 2026-05-28
question: "Are we improving, or just going back and forth?"
verdict: "Surface converging (with one thrash zone); scientific foundation untouched."
method: git log analysis (git rev-list --count = 58)
tags: [meta, audit, revision-history, manuscript]
---

# Revision Trajectory Assessment

**Question:** Over the last ~58 commits, are we getting better, or oscillating?

**One-line verdict:** The surface is mostly *converging*, with one genuine *thrash zone* (Figure 1); but **every one of the 58 commits is packaging-layer — the scientific foundation flagged by the peer-review panel has never been touched.** That is the real source of the "going in circles" feeling.

## Layer 1 — Genuinely converging (monotonic, no reversals)

| Dimension | Evidence | Pattern |
|-----------|----------|---------|
| Word count | 8009 → 6363 → 5006 (v2.5 → v2.6 → v2.7) | Monotonic down to NC cap; never bounced back |
| Structure | Results→Methods offload (`4a24581`, `af6f8fc`) | Done once, stayed |
| Terminology | geometry→topology unification (`9b00837`, `3b1e70b`) | One-time, stable |
| Metadata / reproducibility | ORCID, CITATION.cff, CI auto-release, Zenodo decisions | Additive only — accumulates, never reverts |

This layer shows real, directional progress.

## Layer 2 — Genuine oscillation (the "反反覆覆" feeling)

**Figure 1 alone consumed ~17 of 58 commits**, with explicit revert loops:

- `bb1b0e1` barcode → filtration pedagogy → `acc0aec` biology-first replaces topology pedagogy → `3ac8560` **reverts panel c back to filtration**
- AI 3D render (`00c6ed0`, `77185de`) → later concept-only
- Layout churn: 2-row → full-width → 16:9 → left-align → …

**Root cause:** pixels were adjusted before the figure's *job* ("what question must this panel answer?") was defined. Define the job first; the thrash stops.

## Layer 3 — The foundation was never touched (the decisive finding)

Classifying all 58 commits: **100% are packaging** — word count, figures, ORCID, terminology, submission scaffold. The five gating issues from the 4-reviewer panel (see `peer_review_panel_2026-06-10.md`) have **zero commits** against them:

1. Static loop ≠ temporal cycling (construct validity)
2. max H₁ fragile / no confidence intervals
3. Near-vacuous null; n=50 < own n≥500 standard
4. No negative control
5. YU-006 baseline reported as four different values

From 2026-04-14 to today, not one commit addressed any of these.

## Conclusion

> Not "getting worse," not pure oscillation — **repeatedly repainting a house with a cracked foundation.** The paint (word count, layout, metadata) is genuinely improving. The "going in circles" feeling has two real sources: (a) Figure 1 churned without a defined job; (b) more deeply, the intuition that something load-bearing is unresolved — so no amount of surface editing delivers the feeling of "arrival."

**Break the loop:** the next commit should be the *first one that touches the foundation* — e.g. use the GSE150949 Watermelon lineage barcodes to test whether the H₁ loop is a genuine temporal cycle, resolving reviewer issues #1, #3, and #4 at once. That is worth more than the next ten Figure 1 tweaks.
