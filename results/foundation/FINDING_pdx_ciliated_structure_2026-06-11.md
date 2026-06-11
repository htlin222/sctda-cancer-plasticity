---
title: "The genuine in-vivo topological structure is a drug-emergent ciliated differentiation program"
date: 2026-06-11
analysis: scripts/26_pdx_structure_identity.py + condition/model enrichment
dataset: GSE243562 PDX YU-006 / YU-003, osimertinib residual
status: real biological program; drug-emergent in YU-006; model-specific; origin (transdiff vs expansion) unconfirmed
tags: [positive-finding, pdx, ciliated, transdifferentiation, drug-tolerance]
---

# The PDX in-vivo topological structure is a ciliated differentiation program

Following the controls (the PDX residual carries genuine excess over the dispersion null that
survives cell-cycle removal: G1-only ratio 2.07, 100% bootstrap), we asked what that structure
biologically *is*.

## Result

The structure-driving cells (leave-one-cluster-out: cluster 6, n=89, drops max-H1 by 0.95) are a
coherent **motile-cilia / multiciliated airway-epithelial program**. Top DE genes (vs rest):
PIFO, RSPH1, CFAP126/45/157, CAPS/CAPSL, CETN2, TPPP3, WDR38, LRRC23, FAM81B, C20orf85, C1orf194,
C9orf116, AGR2/AGR3 — an unambiguous ciliogenesis/ciliated signature (with secretory AGR2/AGR3).

## Drug-emergent and model-specific

Scoring the ciliated signature across all four PDX samples (% cells with score > 0.5):

| | untreated | residual |
|---|---:|---:|
| YU-003 | 0.0% | 0.0% |
| YU-006 | 0.1% | **9.6%** |

In YU-006 the program is drug-emergent: in absolute terms ciliated cells rise from ~6 (5,967
cells × 0.1%) to ~187 (1,951 × 9.6%) — a ~30× increase in number even as the tumour compartment
shrinks ~3× under drug. This argues against simple relative enrichment of fixed normal
contaminants (which would keep absolute number roughly constant) and for genuine drug-induced
emergence/expansion. YU-003 shows no such program in either condition — the effect is
**model-specific**.

## Interpretation (appropriately bounded)

EGFR-TKI drug-tolerant persisters are known to adopt differentiated airway/alveolar lineage
states; a drug-emergent multiciliated program in residual disease is consistent with
transdifferentiation toward a tolerant differentiated state. This converts the framework's
"genuine non-cell-cycle in-vivo structure" from an abstract caveat into a concrete, recognisable
biological program.

**Caveats (do not overclaim):**
- Model-specific (1 of 2 PDX). Not a universal drug-tolerance program.
- Definitive proof the ciliated cells are tumour-derived (transdifferentiation) vs an expanded
  normal/non-malignant population requires genotyping/CNV inference not done here. The absolute-
  count argument supports drug-emergence but not cell-of-origin.
- Still fails the traversal test (R=0.965): the population forms a topological feature but does not
  traverse a cycle, so "cyclic plasticity" remains unsupported.

## Bearing on the paper

Strengthens the cautionary/methods paper into a positive one: the control framework does not just
reject confounded signals — it **isolates genuine structure and localises it to a real,
interpretable program**. A focused positive paper ("a drug-emergent ciliated/differentiated
topological state in EGFR-mutant residual disease, validated by dispersion- and cell-cycle-
controlled topology") is viable if the cell-of-origin and a second supporting model can be added.
