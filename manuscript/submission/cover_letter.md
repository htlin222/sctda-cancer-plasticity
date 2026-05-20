# Cover Letter — *Nature Communications* submission

**To:** The Editors, *Nature Communications*
**Subject:** Manuscript submission — "Drug tolerance takes shape: persistent homology reveals cyclic cell-state plasticity in EGFR-mutant lung cancer"

---

Dear Editors,

We are pleased to submit our manuscript, **"Drug tolerance takes shape: persistent homology reveals cyclic cell-state plasticity in EGFR-mutant lung cancer,"** for consideration as a research article in *Nature Communications*.

**The question.** Drug-tolerant persister (DTP) cells survive targeted therapy not by acquiring mutations but by reversibly cycling among transcriptional states. Standard single-cell trajectory inference — Monocle, Slingshot, PAGA, RNA velocity — assumes cells move along tree-like or directed manifolds. That assumption is wrong for reversibly cycling populations and discards the very signal that matters.

**What we did.** We built an open-source pipeline that applies persistent homology (the maximum *H*₁ statistic) and the Mapper algorithm to longitudinal scRNA-seq, with integrated permutation-null testing, two independent cell-cycle controls (Tirosh-gene ablation and stringent S/G2M-score regression), Wasserstein diagram comparison, and Mapper-based gene attribution. We applied it without modification to **five independent EGFR-mutant NSCLC datasets** spanning two TKIs (erlotinib, osimertinib), three biological scales (cell line, PDX, patient tumor), and one 14-patient longitudinal cohort.

**Headline findings.**
- Max *H*₁ rises **monotonically** under osimertinib in three independent systems — PC9 (2.7×), PDX YU-006 (2.2×), PDX YU-003 (1.4×); 4/4 permutation tests pass *p* < 0.01.
- The signal **survives** both Tirosh ablation (ratio 0.71–1.11 across five cohorts) and stringent S/G2M regression (ratio 1.07–4.05, every value > 1); and **survives** Harmony batch correction.
- A 14-patient *Maynard et al.* 2020 EGFR-mutant NSCLC cohort shows monotonic growth across treatment stages (TN 5.37 → PER 7.13 → PD 8.05); 3/3 patients with matched pre/post biopsies increased post-treatment.
- Loop-associated transcripts are enriched for **EMT programs** in both cell line and PDX, with statistically significant but system-specific gene-level overlap — consistent with convergent topology realised through partially different gene repertoires.

**Why *Nature Communications*.** The work sits at the intersection of cancer-cell biology, single-cell genomics, and applied topology — a readership *Nature Communications* serves uniquely well. Methodologically, it introduces the first reproducible, controlled framework for applying TDA to perturbation-response scRNA-seq; biologically, it offers a coordinate-invariant, cross-system phenotype for cell-state plasticity that travels from cell line to PDX to patient tissue. All code is MIT-licensed and the full pipeline is reproducible end-to-end from raw GEO accessions via a single `make pipeline` target.

**Suggested reviewers (TDA-on-biology):**
- Gunnar Carlsson (Stanford) — foundational scTDA methodology
- Raul Rabadan (Columbia) — topological methods for cancer genomics
- Pablo Camara (UPenn) — TDA on single-cell data
- Monica Nicolau (formerly Stanford) — topology of cancer transcriptomes
- Pek Lum (Capella Bio) — Mapper applied to cancer biology

**Suggested reviewers (scRNA-seq + drug tolerance):**
- Aditya Bhatt (Memorial Sloan Kettering) — DTP biology
- Sridhar Ramaswamy / Aaron Hata (MGH) — EGFR drug tolerance
- Yaara Oren (Weizmann) — Watermelon lineage tracing in PC9
- Sydney Shaffer (UPenn) — non-genetic resistance, scRNA-seq

**Excluded reviewers:** none.

**Author contributions** (CRediT): see `credit_taxonomy.md`.
**Competing interests, ethics, data and code availability:** see `declarations.md`.
**Submission checklist:** see `checklist.md`.

This manuscript is not under consideration elsewhere. The accompanying software is publicly archived as GitHub release `v2.7.1-submission` (<https://github.com/htlin222/sctda-cancer-plasticity/releases/tag/v2.7.1-submission>) under the MIT license; the complete analysis pipeline is reproducible end-to-end from raw GEO data via `make pipeline`.

Thank you for considering the manuscript. We look forward to your response.

Sincerely,

**Yueh-Hua Tu, Ph.D.** (corresponding author)
Center for Computational and Systems Biology
National Taiwan University, Taipei, Taiwan
*Email:* yuehhuatu@ntu.edu.tw

on behalf of:

**Hsieh-Ting Lin, M.D.** (first author)
Department of Medical Oncology
Koo Foundation Sun Yat-Sen Cancer Center, Taipei, Taiwan
