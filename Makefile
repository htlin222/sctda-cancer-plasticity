.PHONY: env download preprocess baseline topology ablation biology validate figures pipeline test lint clean pdx kim ablation_perm pdx_genes harmony maynard maynard_full

SHELL := /bin/bash
PYTHON := python
SEED := 42

# ──────────────────────────────────────────────
# Environment
# ──────────────────────────────────────────────
env:
	mamba env create -f envs/environment.yml || conda env create -f envs/environment.yml
	@echo "✅ Run: conda activate sctda"

# ──────────────────────────────────────────────
# Data
# ──────────────────────────────────────────────
download:
	$(PYTHON) scripts/download_data.py

# ──────────────────────────────────────────────
# Analysis Pipeline
# ──────────────────────────────────────────────
preprocess:
	$(PYTHON) scripts/01_preprocess.py --seed $(SEED)

baseline:
	$(PYTHON) scripts/02_standard_analysis.py --seed $(SEED)

topology:
	$(PYTHON) scripts/03_topology.py --seed $(SEED)

ablation:
	$(PYTHON) scripts/04_cell_cycle_ablation.py --seed $(SEED)

biology:
	$(PYTHON) scripts/05_gene_attribution.py --seed $(SEED)

validate:
	$(PYTHON) scripts/06_validation.py --seed $(SEED)

figures:
	$(PYTHON) scripts/07_figures.py

pdx:
	$(PYTHON) scripts/08_validation_pdx.py --seed $(SEED)

kim:
	$(PYTHON) scripts/09_validation_kim_atlas.py --seed $(SEED)

ablation_perm:
	$(PYTHON) scripts/10_ablation_permtest.py --seed $(SEED)

pdx_genes:
	$(PYTHON) scripts/11_pdx_gene_attribution.py --seed $(SEED)

harmony:
	$(PYTHON) scripts/13_harmony_batch_correction.py --seed $(SEED)

maynard:
	$(PYTHON) scripts/14_validation_maynard.py --seed $(SEED)

maynard_full:
	$(PYTHON) scripts/15_maynard_full_cohort.py --seed $(SEED)

# Full pipeline
pipeline: download preprocess baseline topology ablation biology validate pdx kim ablation_perm pdx_genes harmony maynard maynard_full figures
	@echo "✅ Full pipeline complete. Check results/ and manuscript/figures/"

# ──────────────────────────────────────────────
# Quality
# ──────────────────────────────────────────────
test:
	pytest tests/ -v --cov=src/sctda_plasticity

lint:
	ruff check src/ scripts/ tests/
	ruff format --check src/ scripts/ tests/

format:
	ruff format src/ scripts/ tests/

# ──────────────────────────────────────────────
# Cleanup
# ──────────────────────────────────────────────
clean:
	rm -rf data/interim/* data/processed/*
	rm -rf results/figures/* results/tables/*
	@echo "🧹 Cleaned interim/processed data and results"

clean-all: clean
	rm -rf data/raw/*
	@echo "🧹 Cleaned everything including raw data"
