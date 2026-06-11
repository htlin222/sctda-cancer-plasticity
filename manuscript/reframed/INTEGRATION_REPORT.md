# Integration Report — manuscript/reframed/

Consistency edit pass over the seven reframed files, checked against the LOCKED
EVIDENCE table (RESCUE_STRATEGY_2026-06-10.md) and the HARD RULES. Date: 2026-06-10.

Files reviewed: `abstract_title.md`, `introduction.tex`, `results.tex`,
`discussion.tex`, `methods_dispersion_null.tex`, `figure_plan.md`,
`response_to_reviewers.md`.

---

## 1. Surviving cycling / traversal CLAIMS

**None survive as a claim.** Every occurrence of cycling / traversal / reversible
language is framed as (a) the original thesis being withdrawn, (b) a "what we do
not claim" statement, or (c) the circular-coordinate negative control. No file
asserts cells cycle, traverse, or reversibly transition around a loop. Examples
confirming correct (negative) framing:

- `results.tex:90-94` — "we therefore treat this as a built-in negative control
  and explicitly do *not* claim that cells cycle or traverse the loop ... not
  population-wide reversible cycling".
- `discussion.tex:26-27` — "We explicitly do not claim that cells traverse a
  cycle, nor that drug tolerance is a reversible cyclic plasticity."
- `introduction.tex:90-92` — "We therefore describe a non-traversed topological
  reorganisation ... not population-wide cyclic motion."
- `methods_dispersion_null.tex:8,126-128` — explicit "No traversal/cycling claim".
- `abstract_title.md:32` / abstract body — "circular-coordinate analysis confirms
  cells do not traverse it."

**No flags raised in this category.** (One terminology drift, not a claim, is noted
in §3.4.)

---

## 2. Numbers that do not match the locked table

The core max-H1 values, CIs, ratios, gauss-nulls and headline bootstrap fractions
(osi D0/D3/D7/D14; PDX untreated/residual; erl D0/D9/D11; Maynard TN/PER/PD;
R>=0.965 vs 0.881) all match the locked table and are internally arithmetically
consistent (ratio = real/gauss checks out to 2 d.p. throughout). Issues are
confined to **secondary numbers that appear only in some files and are NOT in the
locked table**, plus one internal CV contradiction:

- **2.1 — CV ranges disagree across three files (CONTRADICTION + unsourced).**
  - `methods_dispersion_null.tex:71`: "coefficient of variation 7–19% across the
    osimertinib and PDX groups".
  - `response_to_reviewers.md:125`: osi CVs "D0 5.9%, D3 6.3%, D7 6.8%, D14 8.6%"
    and (126) "PDX residual 23.6%, Maynard PD 44.7%"; (205) "range ~6–45%".
  - `figure_plan.md:67`: "per-cohort CV 6–9%".
  These three ranges are mutually inconsistent (methods floor 7% excludes D0 5.9%
  / D3 6.3%; methods ceiling 19% excludes PDX residual 23.6% and Maynard PD 44.7%;
  figure-plan 6–9% excludes the high-CV groups entirely). None of these CV numbers
  is in the locked evidence table — they are NEW numbers and violate the "no new
  numbers" rule unless traceable to dispersion_null_summary.json.

- **2.2 — Erlotinib bootstrap fractions are new numbers.**
  `response_to_reviewers.md:103,150` report erl D0/D9/D11 "real > Gaussian"
  fractions of **33 / 42 / 0%**. Not in the locked table. (Consistent with the
  locked ratios 0.93/0.99/0.70 in sign, but the % values are unsourced.)

- **2.3 — Maynard bootstrap fractions are new numbers.**
  `response_to_reviewers.md:153` reports Maynard TN/PER/PD "real > Gaussian" of
  **99 / 98 / 95%**. Not in the locked table.

- **2.4 — Rayleigh R = 0.997 (osi D14) is a new number.**
  `figure_plan.md:80-81,93` cite "R = 0.997 (osi D14)". The locked table gives
  only the bound "R >= 0.965 across cohorts". 0.997 is consistent with the bound
  but is a specific new value sourced to circular_coord_summary.json rather than
  the locked table — flag for verification, do not print unless confirmed.

- **2.5 — Watermelon clonal p-value (p = 0.965) is a new number.**
  `figure_plan.md:87,88,95`. Not in the locked table. Also note the numeric
  collision with the Rayleigh bound 0.965 — a reader/typo hazard; verify it is a
  real permutation p-value and not an accidental copy of the R bound.

- **2.6 — Maynard / erlotinib / PDX-untreated gauss-null values shown as "~".**
  `response_to_reviewers.md:150,151,153` and `figure_plan.md:107,111,116-117`
  print approximate nulls ("~1.6", "~3.6", "2.30"). These are back-derivable from
  the locked ratios (erl ≈1.60/1.79/1.54; Maynard ≈3.61; PDX-untr ≈2.31) and are
  consistent, but they are presented as data; ensure they match
  dispersion_null_summary.json exactly before they reach a figure or table.

- **2.7 — S/G2M cell-cycle ablation numbers are new.**
  `figure_plan.md:132-136` (optional Fig 5): "ratio >= 0.71", "PDX YU-006 untreated
  2.89→3.67 (ratio 1.27) and residual 4.65→5.96 (ratio 1.28)". None in the locked
  table. The locked memory note says CC-ablation persistence ratio > 0.7; the 0.71
  is consistent, but the four PDX ablation values are unsourced new numbers.

- **2.8 — Subsample sizes (1,200; n = 379/249) are new.**
  `methods_dispersion_null.tex:77` ("1,200 cells"), `response_to_reviewers.md:184,
  214` ("1,200 ... n = 379/249"). Not in the locked table. Internally consistent
  between methods and response; verify against the analysis script.

---

## 3. Contradictions between sections

- **3.1 — Title mismatch (CONTRADICTION).**
  `abstract_title.md:9` preferred title: "**A dispersion-controlled topological
  signature of drug-tolerant cell-state reorganization in EGFR-mutant lung
  cancer**". `response_to_reviewers.md:4` manuscript_revised title: "**A
  dispersion-controlled topological readout detects drug-induced non-Gaussian
  reorganization of the cell-state space in EGFR-mutant lung cancer**". The two
  files disclose different final titles. Pick one and reconcile.

- **3.2 — PDX null for the residual condition: present in some files, absent in
  methods.** `abstract_title.md:53`, `response_to_reviewers.md:152`,
  `figure_plan.md:108` all give PDX-residual gauss = **2.80** (ratio 1.95 × real
  5.46). `methods_dispersion_null.tex:54-55` states the PDX ratios (1.21, 1.95) but
  never prints the null values, while it does print the osi nulls inline. Minor
  asymmetry; not a contradiction in value, but methods should either print all
  nulls or none for consistency.

- **3.3 — Maynard PER ratio appears in results/discussion/methods but the
  abstract omits patient ratios entirely.** Consistent (abstract intentionally
  scopes out patients), but confirm the abstract's silence is deliberate and that
  the "ratio 1.52–1.69" range in `introduction.tex:96-97` matches the per-stage
  values 1.52/1.69/1.51 used elsewhere (it does; range is correct).

- **3.4 — "12%-occupied control loop" labeled inconsistently (terminology, not
  value).** `figure_plan.md:82` calls it a "synthetic **positive control**";
  `results.tex:86` and `methods_dispersion_null.tex:120` call it "a control loop
  occupied by 12% of cells"; `response_to_reviewers.md:73` says it "genuinely **is
  traversed**"; `introduction.tex:90` says "genuinely **occupied around its
  circumference**". The logic is consistent (the 12% loop is the more-traversed /
  lower-R = 0.881 comparator, real cohorts are more concentrated at R>=0.965), but
  the label oscillates between "control loop", "positive control" and "traversed".
  Standardize to one phrase (suggest: "synthetic comparator loop occupied around
  its circumference (R = 0.881)") to avoid a reviewer reading "positive control"
  as a cycling claim.

- **3.5 — Erlotinib day labels: D9/D11 vs the data table.** All reframed files use
  erl **D0/D9/D11** (matches locked table). Note CLAUDE.md's primary dataset lists
  erlotinib time points D0,D1,D2,D4,D9,D11 — the reframed manuscript reports only
  D0/D9/D11. Confirm this subsetting is intentional and stated in Methods (it is
  not currently explained why D1/D2/D4 are dropped).

- **3.6 — No contradictions found** in: the directional bootstrap fractions
  (D7>D0 100%, D14>D0 100%, PDX residual>untreated 97%) — identical across
  introduction, results, discussion, methods, figure_plan, response. The scope
  statement (osi cell-line + PDX = claim; erlotinib = specificity; patients =
  heterogeneity caveat) is stated consistently in all five prose files.

---

## 4. Punch-list of fixes

1. **Reconcile CV numbers (§2.1).** Choose one authoritative per-group CV set from
   dispersion_null_summary.json. Fix `methods_dispersion_null.tex:71` ("7–19%"),
   `figure_plan.md:67` ("6–9%"), and `response_to_reviewers.md:125-126,205` so all
   three agree. Likely correct global range is ~6–45% (must cover D0 5.9% and
   Maynard PD 44.7%); methods must not say "7–19%".
2. **Source or cut the unsourced secondary numbers (§2.2–2.8):** erl fractions
   33/42/0% · Maynard fractions 99/98/95% · R=0.997 · Watermelon p=0.965 ·
   S/G2M values 2.89→3.67 / 4.65→5.96 / ratios 1.27/1.28/0.71 · subsample sizes
   1,200 and 379/249. Verify each against the JSON/finding files; delete any that
   cannot be traced. (Per the no-new-numbers rule, none of these may be invented.)
3. **Disambiguate p = 0.965 vs Rayleigh R = 0.965 (§2.5).** Confirm the Watermelon
   p-value is genuinely 0.965 and not a transcription of the R bound; reword the
   figure caption so the two 0.965 values cannot be conflated.
4. **Fix the title mismatch (§3.1).** Make `response_to_reviewers.md:4` and the
   preferred title in `abstract_title.md:9` identical.
5. **Standardize the 12%-loop label (§3.4).** One phrase across results,
   methods, intro, figure_plan, response. Drop "positive control" wording in
   `figure_plan.md:82` to avoid implying a demonstrated cycle.
6. **Make methods print the PDX/erl/Maynard nulls or none (§3.2).** Currently osi
   nulls are inline in methods but PDX/erl/Maynard nulls are not; either add them
   (2.31 PDX-untr, 2.80 PDX-res, ~1.60/1.79/1.54 erl, 3.61 Maynard) verified from
   JSON, or move all null values to the single master table and reference it.
7. **State the erlotinib day-subset rationale (§3.5).** Add one sentence in Methods
   explaining why only D0/D9/D11 of the GSE134839 series are analyzed.
8. **Confirm abstract intentionally omits patient ratios (§3.3).** If deliberate
   (scope rule), no change; just verify against final scope statement.

---

### Bottom line
No surviving cycling/traversal claim — the HARD RULES are respected in all prose.
Core locked numbers are clean and self-consistent. All flags are about (i) a real
CV contradiction across three files and (ii) secondary numbers introduced in
`response_to_reviewers.md` and `figure_plan.md` that are not in the locked table
and must be traced to source or removed before submission, plus one title mismatch
and one terminology drift.
