#!/usr/bin/env python3
"""Step 5: Empirical Z-scoring for gene-level AlphaGenome scores."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.append(str(Path(__file__).parent.parent.parent))
from config import PROCESSED_DIR, FIGURES_DIR


def load_universe_file() -> Path | None:
    universe_file = PROCESSED_DIR / "gene_universe.csv"
    if universe_file.exists():
        return universe_file
    print(f"Gene universe file not found: {universe_file}")
    return None


def infer_gene_column(df: pd.DataFrame) -> str | None:
    for col in df.columns:
        low = col.lower()
        if "gene" in low or "symbol" in low:
            return col
    return None


def build_full_gene_matrix(weighted_scores_file: Path, universe_file: Path) -> pd.DataFrame | None:
    if not weighted_scores_file.exists():
        print(f"Weighted score file missing: {weighted_scores_file}")
        return None

    scores_df = pd.read_csv(weighted_scores_file)
    if "gene" not in scores_df.columns:
        print("Weighted score file must include 'gene' column.")
        return None

    uni_df = pd.read_csv(universe_file)
    uni_col = infer_gene_column(uni_df)
    if uni_col is None:
        print("Could not identify gene column in universe file.")
        return None

    score_cols = [c for c in scores_df.columns if c == "weighted_score" or c.startswith("score_")]
    if not score_cols:
        print("No score columns found in weighted score file.")
        return None

    scores_df = scores_df[["gene"] + score_cols].copy()
    for col in score_cols:
        scores_df[col] = pd.to_numeric(scores_df[col], errors="coerce").fillna(0.0)

    all_genes = pd.DataFrame({"gene": pd.Series(uni_df[uni_col].dropna().astype(str).unique())})
    full_df = all_genes.merge(scores_df, on="gene", how="left")
    for col in score_cols:
        full_df[col] = full_df[col].fillna(0.0)

    # Add genes present in scores but not in universe.
    missing_from_universe = scores_df[~scores_df["gene"].isin(set(all_genes["gene"]))]
    if not missing_from_universe.empty:
        full_df = pd.concat([full_df, missing_from_universe], ignore_index=True)

    full_df = full_df.drop_duplicates(subset=["gene"], keep="first")
    print(f"Constructed full gene matrix with {len(full_df)} genes.")
    return full_df


def zscore_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    score_cols = [c for c in out.columns if c == "weighted_score" or c.startswith("score_")]

    for col in score_cols:
        mu = out[col].mean()
        sigma = out[col].std(ddof=0)
        z_col = "z_score" if col == "weighted_score" else f"z_{col.replace('score_', '')}"
        if sigma == 0 or np.isnan(sigma):
            out[z_col] = 0.0
        else:
            out[z_col] = (out[col] - mu) / sigma

    # Keep legacy names for compatibility.
    out["raw_score"] = out["weighted_score"]

    # Empirical P-value based on weighted score ranking (upper tail).
    out = out.sort_values("weighted_score", ascending=False).reset_index(drop=True)
    out["rank"] = out.index + 1
    out["p_empirical"] = out["rank"] / len(out)

    return out


def plot_distribution(scores: pd.Series, output_path: Path) -> None:
    plt.figure(figsize=(10, 6))
    sns.histplot(scores, bins=100, log_scale=(False, True))
    plt.title("Genome-wide Distribution of PP-weighted Regulatory Scores")
    plt.xlabel("Weighted Score")
    plt.ylabel("Count (log scale)")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def main() -> None:
    print("=" * 70)
    print("STEP 5: EMPIRICAL Z-SCORING")
    print("=" * 70)

    scores_file = PROCESSED_DIR / "gene_scores" / "gene_weighted_scores.csv"
    universe_file = load_universe_file()
    if universe_file is None:
        return

    full_df = build_full_gene_matrix(scores_file, universe_file)
    if full_df is None:
        return

    z_df = zscore_columns(full_df)

    out_file = PROCESSED_DIR / "gene_scores" / "gene_z_scores.csv"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    z_df.to_csv(out_file, index=False)
    print(f"Saved calibrated scores to {out_file}")

    fig_path = FIGURES_DIR / "score_distribution_z.png"
    plot_distribution(z_df["weighted_score"], fig_path)
    print(f"Saved distribution plot to {fig_path}")


if __name__ == "__main__":
    main()
