#!/usr/bin/env python3
"""Step 6: Unbiased GSEA using Mann-Whitney U test with FDR correction."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import requests
from scipy import stats
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

sys.path.append(str(Path(__file__).parent.parent.parent))
from config import PROCESSED_DIR, RESULTS_DIR

GO_URL = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2023"
GMT_FILE = PROCESSED_DIR / "GO_Biological_Process_2023.gmt"


def download_gmt_if_needed() -> None:
    if GMT_FILE.exists() and GMT_FILE.stat().st_size > 0:
        print(f"Using cached GMT: {GMT_FILE}")
        return

    print(f"Downloading GMT from {GO_URL}")
    try:
        resp = requests.get(GO_URL, timeout=60)
        resp.raise_for_status()
    except Exception as exc:
        print(f"Failed to download GMT: {exc}")
        return

    GMT_FILE.parent.mkdir(parents=True, exist_ok=True)
    GMT_FILE.write_text(resp.text)
    print(f"Saved GMT to {GMT_FILE}")


def load_gmt(gmt_path: Path, min_genes: int = 15, max_genes: int = 500) -> dict[str, list[str]]:
    pathways: dict[str, list[str]] = {}
    if not gmt_path.exists():
        print(f"GMT not found: {gmt_path}")
        return pathways

    with gmt_path.open("r") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = [g.split(",")[0] for g in parts[2:] if g]
            if min_genes <= len(genes) <= max_genes:
                pathways[name] = genes

    print(f"Loaded {len(pathways)} pathways from {gmt_path.name}")
    return pathways


def run_gsea_for_column(scores_df: pd.DataFrame, score_col: str, pathways: dict[str, list[str]], min_overlap: int = 10) -> pd.DataFrame:
    if "gene" not in scores_df.columns:
        raise ValueError("scores_df must contain 'gene' column")

    work = scores_df[["gene", score_col]].dropna().copy()
    work[score_col] = pd.to_numeric(work[score_col], errors="coerce")
    work = work.dropna(subset=[score_col])

    score_map = dict(zip(work["gene"], work[score_col]))
    available = set(score_map)

    results = []
    for path_name, genes in tqdm(pathways.items(), desc=f"GSEA {score_col}"):
        in_path = [g for g in genes if g in available]
        if len(in_path) < min_overlap:
            continue

        out_path = [g for g in available if g not in set(in_path)]
        if len(out_path) < min_overlap:
            continue

        path_scores = np.array([score_map[g] for g in in_path], dtype=float)
        bg_scores = np.array([score_map[g] for g in out_path], dtype=float)

        if np.all(path_scores == path_scores[0]) and np.all(bg_scores == bg_scores[0]) and path_scores[0] == bg_scores[0]:
            p_val = 1.0
            stat = 0.0
        else:
            stat, p_val = stats.mannwhitneyu(path_scores, bg_scores, alternative="greater")

        results.append(
            {
                "Pathway": path_name,
                "n_genes": len(in_path),
                "p_val": float(p_val),
                "mean_score": float(path_scores.mean()),
                "median_score": float(np.median(path_scores)),
                "stat": float(stat),
            }
        )

    if not results:
        return pd.DataFrame()

    out_df = pd.DataFrame(results)
    out_df["FDR"] = multipletests(out_df["p_val"], method="fdr_bh")[1]
    out_df = out_df.sort_values("p_val").reset_index(drop=True)
    return out_df


def main() -> None:
    print("=" * 70)
    print("STEP 6: MULTI-SCORE GSEA")
    print("=" * 70)

    download_gmt_if_needed()
    pathways = load_gmt(GMT_FILE)
    if not pathways:
        return

    scores_file = PROCESSED_DIR / "gene_scores" / "gene_z_scores.csv"
    if not scores_file.exists():
        print(f"Score file not found: {scores_file}")
        return

    scores_df = pd.read_csv(scores_file)

    score_cols = []
    if "z_score" in scores_df.columns:
        score_cols.append("z_score")
    score_cols.extend(sorted([c for c in scores_df.columns if c.startswith("z_") and c != "z_score"]))

    if not score_cols:
        print("No z-score columns found for enrichment.")
        return

    tables_dir = RESULTS_DIR / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    all_results = []
    for col in score_cols:
        modality = col.replace("z_", "")
        res = run_gsea_for_column(scores_df, col, pathways)
        if res.empty:
            print(f"No pathways passed overlap thresholds for {col}")
            continue

        out_file = tables_dir / f"gsea_enrichment_{modality}.csv"
        res.to_csv(out_file, index=False)
        print(f"Saved {len(res)} rows to {out_file}")

        keep = res.head(50).copy()
        keep["Modality"] = modality
        all_results.append(keep)

    if all_results:
        summary = pd.concat(all_results, ignore_index=True)
        summary_file = tables_dir / "gsea_enrichment_ALL_MODALITIES_summary.csv"
        summary.to_csv(summary_file, index=False)
        print(f"Saved combined summary to {summary_file}")


if __name__ == "__main__":
    main()
