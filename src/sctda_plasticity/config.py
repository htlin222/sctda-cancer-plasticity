"""
Project-wide configuration and constants.
"""

from pathlib import Path

# ──────────────────────────────────────────────
# Paths
# ──────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_INTERIM = PROJECT_ROOT / "data" / "interim"
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
TABLES_DIR = RESULTS_DIR / "tables"

# ──────────────────────────────────────────────
# Analysis parameters
# ──────────────────────────────────────────────
SEED = 42

# Preprocessing
MIN_GENES = 200
MAX_GENES = 5000
MAX_PCT_MITO = 20.0
MIN_CELLS = 3
N_TOP_GENES = 3000
N_PCS = 50

# Persistent homology
PH_MAX_DIM = 1        # H₀ and H₁
PH_N_PCS = 30         # PCs used for PH
N_PERMUTATIONS = 500  # Permutation test

# Mapper
MAPPER_N_CUBES = 15
MAPPER_OVERLAP = 0.3

# ──────────────────────────────────────────────
# Datasets
# ──────────────────────────────────────────────
DATASETS = {
    "discovery": {
        "accession": "GSE134839",
        "name": "PC9 erlotinib time-series",
        "reference": "Aissa et al., Nat Commun 2021",
        "timepoints": ["D0", "D1", "D2", "D4", "D9", "D11"],
        "samples": {
            "GSM3972657": "D0",
            "GSM3972658": "D1",
            "GSM3972659": "D2",
            "GSM3972660": "D4",
            "GSM3972661": "D9",
            "GSM3972662": "D11",
        },
    },
    "validation_1": {
        "accession": "GSE150949",
        "name": "PC9 osimertinib",
        "reference": "GEO",
        "timepoints": None,
    },
}

# ──────────────────────────────────────────────
# Visualization
# ──────────────────────────────────────────────
TIMEPOINT_COLORS = {
    "D0": "#2196F3",
    "D1": "#4CAF50",
    "D2": "#FF9800",
    "D4": "#FF5722",
    "D9": "#F44336",
    "D11": "#9C27B0",
}

FIGURE_DPI = 300
FIGURE_FORMAT = ["pdf", "png"]
