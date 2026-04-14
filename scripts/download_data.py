#!/usr/bin/env python
"""
Download datasets from GEO for the scTDA cancer plasticity project.

GSE134839: 6 Drop-seq DGE files (D0, D1, D2, D4, D9, D11), ~848 cells total
GSE150949: PC9 count matrix + metadata CSV files

Usage:
    python scripts/download_data.py
    python scripts/download_data.py --dataset GSE134839
"""

import argparse
import gzip
import logging
import subprocess
import urllib.request
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

DATA_DIR = Path("data/raw")

# ── GSE134839: Drop-seq DGE files ───────────────────────────
GSE134839_FILES = {
    "GSM3972657_D0.dge.txt.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3972nnn/GSM3972657/suppl/GSM3972657_D0.dge.txt.gz",
        "timepoint": "D0",
    },
    "GSM3972658_D1.dge.txt.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3972nnn/GSM3972658/suppl/GSM3972658_D1.dge.txt.gz",
        "timepoint": "D1",
    },
    "GSM3972659_D2.dge.txt.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3972nnn/GSM3972659/suppl/GSM3972659_D2.dge.txt.gz",
        "timepoint": "D2",
    },
    "GSM3972660_D4.dge.txt.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3972nnn/GSM3972660/suppl/GSM3972660_D4.dge.txt.gz",
        "timepoint": "D4",
    },
    "GSM3972661_D9.dge.txt.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3972nnn/GSM3972661/suppl/GSM3972661_D9.dge.txt.gz",
        "timepoint": "D9",
    },
    "GSM3972662_D11.dge.txt.gz": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3972nnn/GSM3972662/suppl/GSM3972662_D11.dge.txt.gz",
        "timepoint": "D11",
    },
}

# ── GSE150949: PC9 osimertinib ──────────────────────────────
GSE150949_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150949/suppl"
GSE150949_FILES = {
    "GSE150949_pc9_count_matrix.csv.gz": {
        "url": f"{GSE150949_BASE}/GSE150949_pc9_count_matrix.csv.gz",
    },
    "GSE150949_metaData_with_lineage.txt.gz": {
        "url": f"{GSE150949_BASE}/GSE150949_metaData_with_lineage.txt.gz",
    },
}


def download_file(url: str, dest: Path, retries: int = 3) -> bool:
    """Download a file with retry logic."""
    if dest.exists() and dest.stat().st_size > 0:
        logger.info(f"  Already exists: {dest.name}")
        return True

    for attempt in range(retries):
        try:
            logger.info(f"  Downloading: {dest.name} (attempt {attempt + 1})")
            urllib.request.urlretrieve(url, dest)
            logger.info(f"  Downloaded: {dest.name} ({dest.stat().st_size:,} bytes)")
            return True
        except Exception as e:
            logger.warning(f"  Attempt {attempt + 1} failed: {e}")
            if dest.exists():
                dest.unlink()

    # Fallback to curl
    try:
        logger.info(f"  Trying curl for {dest.name}...")
        result = subprocess.run(
            ["curl", "-L", "-o", str(dest), url],
            capture_output=True, text=True, timeout=300,
        )
        if result.returncode == 0 and dest.exists() and dest.stat().st_size > 0:
            logger.info(f"  Downloaded via curl: {dest.name}")
            return True
    except Exception:
        pass

    logger.error(f"  FAILED to download: {dest.name}")
    return False


def verify_dge_file(filepath: Path) -> dict:
    """Verify a Drop-seq DGE file and report stats."""
    try:
        opener = gzip.open if str(filepath).endswith(".gz") else open
        with opener(filepath, "rt") as f:
            header = f.readline().strip().split("\t")
            n_cells = len(header) - 1  # first column is gene name
            n_genes = sum(1 for _ in f)
        return {"n_cells": n_cells, "n_genes": n_genes, "valid": True}
    except Exception as e:
        return {"error": str(e), "valid": False}


def download_gse134839(output_dir: Path) -> None:
    """Download GSE134839 Drop-seq DGE files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    success = 0
    total_cells = 0

    for filename, info in GSE134839_FILES.items():
        dest = output_dir / filename
        if download_file(info["url"], dest):
            stats = verify_dge_file(dest)
            if stats.get("valid"):
                logger.info(
                    f"  Verified {info['timepoint']}: "
                    f"{stats['n_cells']} cells × {stats['n_genes']} genes"
                )
                total_cells += stats["n_cells"]
                success += 1
            else:
                logger.warning(f"  Verification failed for {filename}: {stats.get('error')}")

    logger.info(f"\nGSE134839: {success}/{len(GSE134839_FILES)} files downloaded")
    logger.info(f"Total cells across timepoints: {total_cells}")


def download_gse150949(output_dir: Path) -> None:
    """Download GSE150949 PC9 osimertinib files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    for filename, info in GSE150949_FILES.items():
        dest = output_dir / filename
        download_file(info["url"], dest)


def main():
    parser = argparse.ArgumentParser(description="Download GEO datasets")
    parser.add_argument(
        "--dataset",
        choices=["GSE134839", "GSE150949", "all"],
        default="all",
        help="Which dataset to download",
    )
    args = parser.parse_args()

    DATA_DIR.mkdir(parents=True, exist_ok=True)

    if args.dataset in ("GSE134839", "all"):
        logger.info(f"\n{'='*60}")
        logger.info("Dataset: GSE134839 — PC9 erlotinib time-series")
        logger.info("Reference: Aissa et al., Nat Commun 2021")
        logger.info(f"{'='*60}")
        download_gse134839(DATA_DIR / "GSE134839")

    if args.dataset in ("GSE150949", "all"):
        logger.info(f"\n{'='*60}")
        logger.info("Dataset: GSE150949 — PC9 osimertinib (validation)")
        logger.info(f"{'='*60}")
        download_gse150949(DATA_DIR / "GSE150949")

    logger.info("\nDownload complete. Check data/raw/ for files.")


if __name__ == "__main__":
    main()
