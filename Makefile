.PHONY: env download preprocess baseline topology ablation biology validate figures pipeline test lint clean

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

# Full pipeline
pipeline: download preprocess baseline topology ablation biology validate figures
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
