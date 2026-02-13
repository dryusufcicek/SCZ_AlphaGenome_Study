#!/usr/bin/env python3
"""Generate key manuscript figures and core supplementary tables."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.append(str(Path(__file__).parent.parent))
from config import PROCESSED_DIR, FIGURES_DIR, RESULTS_DIR, CT_DIR

FIGURES_DIR.mkdir(parents=True, exist_ok=True)
TABLES_DIR = RESULTS_DIR / "tables"
TABLES_DIR.mkdir(parents=True, exist_ok=True)


def generate_table_s1_variants() -> None:
    cs_file = PROCESSED_DIR / "credible_sets" / "scz_credible_sets_official.csv"
    if not cs_file.exists():
        print("Skipping Table S1: credible set file missing")
        return
    pd.read_csv(cs_file).to_csv(TABLES_DIR / "Table_S1_Variants.csv", index=False)


def generate_table_s2_genes() -> pd.DataFrame | None:
    z_file = PROCESSED_DIR / "gene_scores" / "gene_z_scores.csv"
    if not z_file.exists():
        print("Skipping Table S2: z-score file missing")
        return None
    df = pd.read_csv(z_file)
    df.to_csv(TABLES_DIR / "Table_S2_Gene_Analysis.csv", index=False)
    return df


def generate_table_s4_celltype() -> pd.DataFrame | None:
    ct_file = CT_DIR / "cell_type_footprint_enrichment.csv"
    if not ct_file.exists():
        print("Skipping Table S4: cell-type file missing")
        return None
    df = pd.read_csv(ct_file)
    df.to_csv(TABLES_DIR / "Table_S4_CellType_Specificity.csv", index=False)
    return df


def plot_fig2a_regulatory_landscape(gene_df: pd.DataFrame) -> None:
    if "z_score" not in gene_df.columns:
        print("Skipping Fig2A: z_score column missing")
        return

    plot_df = gene_df.sort_values("z_score", ascending=False).reset_index(drop=True)
    top_genes = plot_df.head(5)

    plt.figure(figsize=(10, 6))
    plt.plot(plot_df.index, plot_df["z_score"], color="#1f3b73", linewidth=1.5)
    plt.scatter(top_genes.index, top_genes["z_score"], color="#c0392b", zorder=3)

    for i, row in top_genes.iterrows():
        plt.text(i + max(10, len(plot_df) // 100), row["z_score"], str(row["gene"]), fontsize=9)

    plt.title("Genome-wide Regulatory Burden Ranking")
    plt.xlabel("Gene Rank")
    plt.ylabel("Empirical Z-score")
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "Fig2A_Regulatory_Landscape.png", dpi=300)
    plt.close()


def plot_fig2c_sensitivity() -> None:
    primary_file = TABLES_DIR / "gsea_enrichment_score.csv"
    sens_file = TABLES_DIR / "sensitivity_results_no_top1pct_no_mhc.csv"
    if not primary_file.exists() or not sens_file.exists():
        print("Skipping Fig2C: sensitivity input files missing")
        return

    primary = pd.read_csv(primary_file).head(10).copy()
    sens = pd.read_csv(sens_file).head(10).copy()

    for df in (primary, sens):
        if "p_val" not in df.columns:
            p_alias = next((c for c in ["P_Value", "P-value"] if c in df.columns), None)
            if p_alias is not None:
                df["p_val"] = df[p_alias]

    if "p_val" not in primary.columns or "p_val" not in sens.columns:
        print("Skipping Fig2C: p-value column missing")
        return

    primary["neglog10p"] = -np.log10(primary["p_val"].clip(lower=1e-300))
    sens["neglog10p"] = -np.log10(sens["p_val"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(data=primary, x="neglog10p", y="Pathway", color="#4f81bd", alpha=0.75, ax=ax)
    sns.barplot(data=sens, x="neglog10p", y="Pathway", color="#e74c3c", alpha=0.45, ax=ax)
    ax.set_title("Pathway Signal Robustness Under Sensitivity Filters")
    ax.set_xlabel("-log10(P)")
    ax.set_ylabel("Pathway")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "Fig2C_Sensitivity_Comparison.png", dpi=300)
    plt.close()


def plot_fig4c_celltype(ct_df: pd.DataFrame | None) -> None:
    if ct_df is None or ct_df.empty:
        print("Skipping Fig4C: no cell-type dataframe")
        return

    work = ct_df.rename(columns={"P_Binom": "p_val", "Cell": "label", "CellType": "label", "P_Norm": "p_val"}).copy()
    if "p_val" not in work.columns or "label" not in work.columns:
        print("Skipping Fig4C: required columns missing")
        return

    work["neglog10p"] = -np.log10(pd.to_numeric(work["p_val"], errors="coerce").clip(lower=1e-300))

    plt.figure(figsize=(9, 6))
    sns.barplot(data=work, x="label", y="neglog10p", color="#2e8b57")
    plt.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=1)
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("-log10(P)")
    plt.title("Cell-type Specificity (Footprint-normalized)")
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "Fig4C_CellType_Specificity.png", dpi=300)
    plt.close()


def main() -> None:
    print("Generating manuscript tables and figures...")
    generate_table_s1_variants()
    gene_df = generate_table_s2_genes()
    ct_df = generate_table_s4_celltype()

    if gene_df is not None:
        plot_fig2a_regulatory_landscape(gene_df)
    plot_fig2c_sensitivity()
    plot_fig4c_celltype(ct_df)


if __name__ == "__main__":
    main()
