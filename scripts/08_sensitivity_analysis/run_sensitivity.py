#!/usr/bin/env python3
"""Step 8: Sensitivity analysis for pathway enrichment robustness."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pandas as pd

sys.path.append(str(Path(__file__).parent.parent.parent))
from config import PROCESSED_DIR, RESULTS_DIR

# Reuse Step 6 implementation.
gsea_script = Path(__file__).parent.parent / "06_enrichment" / "run_gsea.py"
spec = importlib.util.spec_from_file_location("gsea_module", gsea_script)
gsea_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gsea_module)


def _exclude_mhc_genes(df: pd.DataFrame) -> pd.DataFrame:
    mhc_pattern = r"^HLA-|^C4[AB]|^MICA|^MICB|^TNF|^LTA"
    mask = df["gene"].astype(str).str.contains(mhc_pattern, na=False, regex=True)
    removed = int(mask.sum())
    print(f"Excluded {removed} genes by MHC pattern filter.")
    return df.loc[~mask].copy()


def run_sensitivity_analysis(exclude_top_1pct: bool, exclude_mhc: bool) -> None:
    print("=" * 70)
    print(f"Sensitivity run: exclude_top_1pct={exclude_top_1pct}, exclude_mhc={exclude_mhc}")
    print("=" * 70)

    z_file = PROCESSED_DIR / "gene_scores" / "gene_z_scores.csv"
    if not z_file.exists():
        print(f"Missing Z-score file: {z_file}")
        return

    df = pd.read_csv(z_file)
    if "gene" not in df.columns or "z_score" not in df.columns:
        print("Z-score file must contain gene and z_score columns.")
        return

    work = df[["gene", "z_score"]].copy()

    if exclude_top_1pct:
        cutoff = work["z_score"].quantile(0.99)
        work = work[work["z_score"] <= cutoff].copy()
        print(f"Excluded top 1% by z_score cutoff {cutoff:.4f}")

    if exclude_mhc:
        work = _exclude_mhc_genes(work)

    pathways = gsea_module.load_gmt(gsea_module.GMT_FILE)
    if not pathways:
        gsea_module.download_gmt_if_needed()
        pathways = gsea_module.load_gmt(gsea_module.GMT_FILE)
    if not pathways:
        print("No pathways available.")
        return

    res = gsea_module.run_gsea_for_column(work, "z_score", pathways)
    if res.empty:
        print("No pathway results generated.")
        return

    suffix = ""
    if exclude_top_1pct:
        suffix += "_no_top1pct"
    if exclude_mhc:
        suffix += "_no_mhc"

    out_dir = RESULTS_DIR / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"sensitivity_results{suffix}.csv"
    res.to_csv(out_file, index=False)
    print(f"Saved sensitivity results to {out_file}")


def main() -> None:
    run_sensitivity_analysis(exclude_top_1pct=True, exclude_mhc=False)
    run_sensitivity_analysis(exclude_top_1pct=False, exclude_mhc=True)
    run_sensitivity_analysis(exclude_top_1pct=True, exclude_mhc=True)


if __name__ == "__main__":
    main()
