#!/usr/bin/env python3
"""Step 9a: Hi-C chromatin loop overlap validation."""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

sys.path.append(str(Path(__file__).parent.parent.parent))
from config import PROCESSED_DIR, RESULTS_DIR

BASE_DIR = Path(__file__).resolve().parents[2]

HIC_CANDIDATES = [
    BASE_DIR / "data" / "validation" / "adultbrain_hic.bedpe",
    BASE_DIR.parent / "data" / "validation" / "adultbrain_hic.bedpe",
]


def resolve_hic_file() -> Path | None:
    for candidate in HIC_CANDIDATES:
        if candidate.exists():
            return candidate
    return None


def load_snps() -> pd.DataFrame | None:
    snps_file = PROCESSED_DIR / "credible_sets" / "scz_credible_sets_official.csv"
    if not snps_file.exists():
        print(f"Credible sets not found: {snps_file}")
        return None
    df = pd.read_csv(snps_file)
    return df[["SNP", "CHR", "BP"]].copy()


def load_hic_loops(hic_file: Path) -> pd.DataFrame | None:
    try:
        df = pd.read_csv(hic_file, sep=r"\s+")
    except Exception as exc:
        print(f"Could not read Hi-C file: {exc}")
        return None

    if "chrom1" not in df.columns and len(df.columns) >= 6:
        df.columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"] + list(df.columns[6:])

    needed = {"chrom1", "start1", "end1", "chrom2", "start2", "end2"}
    if not needed.issubset(df.columns):
        print(f"Hi-C file missing required columns: {needed - set(df.columns)}")
        return None

    return df


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ms, me = merged[-1]
        if s <= me:
            merged[-1] = (ms, max(me, e))
        else:
            merged.append((s, e))
    return merged


def build_anchor_index(loops_df: pd.DataFrame) -> dict[str, list[tuple[int, int]]]:
    anchors: dict[str, list[tuple[int, int]]] = {}
    for _, row in loops_df.iterrows():
        for chrom_col, start_col, end_col in [("chrom1", "start1", "end1"), ("chrom2", "start2", "end2")]:
            chrom = str(row[chrom_col])
            anchors.setdefault(chrom, []).append((int(row[start_col]), int(row[end_col])))

    for chrom in list(anchors):
        anchors[chrom] = merge_intervals(anchors[chrom])
    return anchors


def snp_hits_anchor(chrom: str, pos: int, anchors: dict[str, list[tuple[int, int]]]) -> bool:
    intervals = anchors.get(chrom)
    if not intervals:
        return False
    for start, end in intervals:
        if start <= pos <= end:
            return True
        if start > pos:
            break
    return False


def main() -> None:
    print("=" * 70)
    print("STEP 9a: Hi-C LOOP OVERLAP VALIDATION")
    print("=" * 70)

    hic_file = resolve_hic_file()
    if hic_file is None:
        print("No Hi-C BEDPE file found in expected locations.")
        return

    snps_df = load_snps()
    if snps_df is None:
        return

    loops_df = load_hic_loops(hic_file)
    if loops_df is None:
        return

    anchors = build_anchor_index(loops_df)
    print(f"Testing {len(snps_df)} SNPs against {len(loops_df)} loops ({hic_file})")

    rows = []
    hits = 0
    for _, row in snps_df.iterrows():
        chrom = str(row["CHR"])
        chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
        pos = int(row["BP"])
        in_loop = snp_hits_anchor(chrom, pos, anchors)
        if in_loop:
            hits += 1
        rows.append({"SNP": row["SNP"], "CHR": row["CHR"], "BP": pos, "InLoop": in_loop})

    rate = 100.0 * hits / len(snps_df)

    total_anchor_bp = 0
    for intervals in anchors.values():
        total_anchor_bp += sum(e - s for s, e in intervals)

    expected_rate = 100.0 * total_anchor_bp / 3.1e9
    enrichment = rate / expected_rate if expected_rate > 0 else float("inf")

    print(f"Observed in-loop SNPs: {hits}/{len(snps_df)} ({rate:.2f}%)")
    print(f"Expected anchor coverage: {expected_rate:.2f}%")
    print(f"Enrichment: {enrichment:.2f}x")

    out_dir = RESULTS_DIR / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "hic_loop_overlap.csv"
    pd.DataFrame(rows).query("InLoop == True").to_csv(out_file, index=False)
    print(f"Saved hit table to {out_file}")


if __name__ == "__main__":
    main()
