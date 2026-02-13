#!/usr/bin/env python3
"""Repair and normalize AlphaGenome variant score CSV schema.

This script preserves all rows from mixed historical formats and writes a
canonical 12-column file expected by downstream scripts.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import List, Optional

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data" / "processed"
IN_FILE = DATA_DIR / "alphagenome_variant_scores.csv"
BACKUP_FILE = DATA_DIR / "alphagenome_variant_scores.csv.bak_schema_fix"
CREDIBLE_SET_FILE = DATA_DIR / "credible_sets" / "scz_credible_sets_official.csv"

CANONICAL_HEADER = [
    "SNP",
    "target_gene",
    "target_gene_score",
    "all_genes",
    "success",
    "score_DNase",
    "score_H3K27ac",
    "score_H3K4me1",
    "score_H3K4me3",
    "score_CAGE",
    "score_CTCF",
    "score_RNA",
]


def _is_success_token(value: str) -> bool:
    return value.strip().lower() in {"true", "false", "1", "0"}


def normalize_row(row: List[str]) -> Optional[List[str]]:
    """Map variable-width historical row into canonical 12-column schema."""
    if not row:
        return None

    if row[0] == "SNP":
        return None

    if len(row) < 5:
        return None

    snp = row[0]
    target_gene = row[1] if len(row) > 1 else ""
    target_gene_score = row[2] if len(row) > 2 else ""

    success_idx = None
    for i in range(len(row) - 1, 2, -1):
        if _is_success_token(row[i]):
            success_idx = i
            break

    if success_idx is None:
        # Fallback for unexpected rows: assume legacy 5-column layout.
        success_idx = 4 if len(row) >= 5 else len(row) - 1

    all_genes = ",".join(row[3:success_idx]) if success_idx > 3 else ""
    success = row[success_idx] if success_idx < len(row) else ""

    modalities = row[success_idx + 1 :]
    if len(modalities) > 7:
        modalities = modalities[:7]
    if len(modalities) < 7:
        modalities += [""] * (7 - len(modalities))

    out = [snp, target_gene, target_gene_score, all_genes, success] + modalities
    if len(out) != 12:
        return None
    return out


def main() -> None:
    print("=" * 72)
    print("FIXING CSV SCHEMA: alphagenome_variant_scores.csv")
    print("=" * 72)

    if not IN_FILE.exists() and not BACKUP_FILE.exists():
        print(f"Error: neither {IN_FILE} nor {BACKUP_FILE} exists.")
        return

    source_file = BACKUP_FILE if BACKUP_FILE.exists() else IN_FILE

    # Keep a stable backup once.
    if IN_FILE.exists() and not BACKUP_FILE.exists():
        print(f"Creating backup: {BACKUP_FILE}")
        IN_FILE.rename(BACKUP_FILE)
        source_file = BACKUP_FILE

    print(f"Reading source rows from: {source_file}")

    normalized = []
    dropped = 0
    seen = set()

    with open(source_file, "r", newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            fixed = normalize_row(row)
            if fixed is None:
                dropped += 1
                continue

            snp = fixed[0]
            if snp in seen:
                # Keep first occurrence to avoid accidental duplication across reruns.
                continue
            seen.add(snp)
            normalized.append(fixed)

    print(f"Recovered rows: {len(normalized)}")
    print(f"Dropped unparsable/header rows: {dropped}")

    with open(IN_FILE, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(CANONICAL_HEADER)
        writer.writerows(normalized)

    print(f"Wrote canonical file: {IN_FILE}")

    if CREDIBLE_SET_FILE.exists():
        import pandas as pd

        cs = pd.read_csv(CREDIBLE_SET_FILE, usecols=["SNP"])
        out = pd.read_csv(IN_FILE, usecols=["SNP"])
        missing = set(cs["SNP"]) - set(out["SNP"])
        print(f"Coverage vs official credible set: {len(out)} / {len(cs)}")
        print(f"Missing SNPs after repair: {len(missing)}")


if __name__ == "__main__":
    main()
