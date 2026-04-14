#!/usr/bin/env python
"""
Convert large CSV tables to Parquet for git repo storage efficiency.

Keeps CSV for small files (<50KB) for human readability; converts larger
gene/DE tables to parquet (typically 3-10x smaller).

Run: python scripts/12_convert_tables_to_parquet.py
"""

import logging
import sys
from pathlib import Path

sys.path.insert(0, "src")

import pandas as pd

from sctda_plasticity.config import TABLES_DIR

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# Threshold: convert CSVs larger than this to parquet
SIZE_THRESHOLD_BYTES = 50_000


def main():
    total_csv_bytes = 0
    total_parquet_bytes = 0
    converted = []

    for csv_path in sorted(TABLES_DIR.glob("*.csv")):
        size = csv_path.stat().st_size
        total_csv_bytes += size

        if size < SIZE_THRESHOLD_BYTES:
            logger.info(f"  SKIP (small): {csv_path.name} ({size:,} B)")
            continue

        parquet_path = csv_path.with_suffix(".parquet")
        try:
            df = pd.read_csv(csv_path)
            df.to_parquet(parquet_path, compression="snappy", index=False)
            p_size = parquet_path.stat().st_size
            ratio = p_size / size
            total_parquet_bytes += p_size
            converted.append((csv_path.name, size, p_size, ratio))
            logger.info(
                f"  CONVERT: {csv_path.name}: {size:,} B -> "
                f"{parquet_path.name}: {p_size:,} B ({ratio:.1%})"
            )
        except Exception as e:
            logger.warning(f"  FAILED: {csv_path.name}: {e}")

    logger.info("\n=== Summary ===")
    logger.info(f"Converted {len(converted)} file(s)")
    if converted:
        saved = sum(c[1] - c[2] for c in converted)
        logger.info(f"Original size: {sum(c[1] for c in converted):,} B")
        logger.info(f"Parquet size:  {sum(c[2] for c in converted):,} B")
        logger.info(f"Saved:         {saved:,} B "
                    f"({saved / sum(c[1] for c in converted):.1%})")

    logger.info("\nNote: CSV files retained for readability; parquet files added.")
    logger.info("Add CSVs > 50KB to .gitignore if you want git to track only parquet.")


if __name__ == "__main__":
    main()
