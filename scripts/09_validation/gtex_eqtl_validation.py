#!/usr/bin/env python3
"""Step 9b: GTEx eQTL cross-validation for AlphaGenome predictions.

This script performs a best-effort validation. If GTEx files are absent,
it writes an explicit stub report rather than silently failing.
"""

from __future__ import annotations

from pathlib import Path
import gzip

import pandas as pd
import requests

import sys

sys.path.append(str(Path(__file__).parent.parent.parent))
from config import PROCESSED_DIR, RESULTS_DIR

BASE_DIR = Path(__file__).resolve().parents[2]
GTEX_DIR = BASE_DIR / "data" / "validation" / "gtex"
GTEX_DIR.mkdir(parents=True, exist_ok=True)

TARGET_TISSUES = [
    "Brain_Frontal_Cortex_BA9",
    "Brain_Cortex",
]


def download_gtex_brain_eqtls() -> list[Path]:
    base_url = "https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL/"
    downloaded = []

    for tissue in TARGET_TISSUES:
        filename = f"{tissue}.v8.signif_variant_gene_pairs.txt.gz"
        out = GTEX_DIR / filename

        if out.exists() and out.stat().st_size > 1024:
            downloaded.append(out)
            continue

        url = base_url + filename
        print(f"Downloading {tissue} GTEx file...")
        try:
            resp = requests.get(url, stream=True, timeout=60)
            if resp.status_code != 200:
                print(f"Failed download for {tissue}: HTTP {resp.status_code}")
                continue
            with open(out, "wb") as f:
                for chunk in resp.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
            downloaded.append(out)
            print(f"Saved {out}")
        except Exception as exc:
            print(f"Download error for {tissue}: {exc}")

    return downloaded


def load_gtex_eqtl_file(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    try:
        with gzip.open(path, "rt") as f:
            first = f.readline()
            if first.startswith("<?xml"):
                return None
        return pd.read_csv(path, sep="\t", compression="gzip")
    except Exception:
        return None


def build_stub_report(ag_scores: pd.DataFrame, reason: str) -> pd.DataFrame:
    high_effect_n = 0
    if "target_gene_score" in ag_scores.columns:
        threshold = ag_scores["target_gene_score"].abs().quantile(0.75)
        high_effect_n = int((ag_scores["target_gene_score"].abs() > threshold).sum())

    return pd.DataFrame(
        [
            {
                "analysis": "GTEx eQTL validation",
                "status": "NOT_COMPLETED",
                "reason": reason,
                "alphagenome_variants": int(len(ag_scores)),
                "high_effect_variants": high_effect_n,
                "variants_validated_by_eqtl": pd.NA,
                "validation_rate": pd.NA,
            }
        ]
    )


def main() -> None:
    print("=" * 70)
    print("STEP 9b: GTEx eQTL VALIDATION")
    print("=" * 70)

    ag_file = PROCESSED_DIR / "alphagenome_variant_scores.csv"
    if not ag_file.exists():
        print(f"Missing AlphaGenome score file: {ag_file}")
        return

    ag_scores = pd.read_csv(ag_file)
    print(f"Loaded {len(ag_scores)} AlphaGenome variant rows")

    downloaded = download_gtex_brain_eqtls()
    eqtl_df = None
    for path in downloaded:
        eqtl_df = load_gtex_eqtl_file(path)
        if eqtl_df is not None and not eqtl_df.empty:
            print(f"Loaded {len(eqtl_df)} rows from {path.name}")
            break

    out_dir = RESULTS_DIR / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)

    if eqtl_df is None:
        stub = build_stub_report(ag_scores, "GTEx file unavailable or unreadable")
        out = out_dir / "gtex_validation_stub.csv"
        stub.to_csv(out, index=False)
        print(f"Saved stub report to {out}")
        return

    # Coordinate-level harmonization is dataset-specific and not completed in this script.
    # We save a transparent placeholder summary instead of reporting unsupported concordance.
    stub = build_stub_report(ag_scores, "Coordinate harmonization not implemented in this version")
    stub.loc[0, "status"] = "PARTIAL"
    out = out_dir / "gtex_validation_stub.csv"
    stub.to_csv(out, index=False)
    print(f"Saved partial validation report to {out}")


if __name__ == "__main__":
    main()
