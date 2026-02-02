#!/usr/bin/env python3
"""
03_pathway_enrichment.py
========================

Step 3: Ranked Gene Set Enrichment Analysis (GSEA) for pathway enrichment.

This script tests whether neurobiologically-relevant pathways are enriched
among high-scoring AlphaGenome targets using a rank-based approach.

Input:
    - data/processed/alphagenome_scores.csv

Output:
    - results/tables/gsea_results.csv
    - results/figures/gsea_enrichment.png

Usage:
    python scripts/03_pathway_enrichment.py

Author: Yusuf Cicek
Date: February 2026
"""

import pandas as pd
import numpy as np
from scipy import stats
from collections import defaultdict
from pathlib import Path
import matplotlib.pyplot as plt
import sys

sys.path.append(str(Path(__file__).parent.parent))
from config import PROCESSED_DIR, TABLES_DIR, FIGURES_DIR, PATHWAY_SETS


def load_ranked_genes(scores_file):
    """
    Create ranked gene list from AlphaGenome scores.
    
    Parameters
    ----------
    scores_file : Path
        Path to AlphaGenome scores CSV
        
    Returns
    -------
    pd.DataFrame
        Genes ranked by effect score (highest first)
    """
    scores = pd.read_csv(scores_file)
    
    gene_scores = defaultdict(list)
    
    for _, row in scores.iterrows():
        if pd.notna(row.get('genes')) and row['genes']:
            for entry in str(row['genes']).split(';'):
                if ':' in entry:
                    gene, score = entry.split(':')
                    gene_scores[gene].append(abs(float(score)))
    
    ranked = pd.DataFrame([
        {'gene': g, 'score': max(s)}
        for g, s in gene_scores.items()
    ]).sort_values('score', ascending=False).reset_index(drop=True)
    
    ranked['rank'] = range(1, len(ranked) + 1)
    
    return ranked


def gsea_mannwhitney(ranked_genes, pathway_genes):
    """
    Perform GSEA using Mann-Whitney U test.
    
    Tests whether pathway genes have lower ranks (higher scores)
    than non-pathway genes.
    """
    pathway_set = set(pathway_genes)
    
    pathway_ranks = ranked_genes[ranked_genes['gene'].isin(pathway_set)]['rank'].values
    non_pathway_ranks = ranked_genes[~ranked_genes['gene'].isin(pathway_set)]['rank'].values
    
    n_found = len(pathway_ranks)
    
    if n_found < 2:
        return None, None, None, 0, n_found
    
    statistic, p_value = stats.mannwhitneyu(
        pathway_ranks, non_pathway_ranks, 
        alternative='less'
    )
    
    # Rank biserial correlation (effect size)
    n1, n2 = len(pathway_ranks), len(non_pathway_ranks)
    r = 1 - (2 * statistic) / (n1 * n2)
    
    median_rank = np.median(pathway_ranks)
    
    return statistic, p_value, r, median_rank, n_found


def main():
    """Main GSEA pipeline."""
    print("=" * 70)
    print("STEP 3: PATHWAY ENRICHMENT ANALYSIS (GSEA)")
    print("=" * 70)
    
    # Load ranked genes
    scores_file = PROCESSED_DIR / "alphagenome_scores.csv"
    if not scores_file.exists():
        print(f"Error: Scores not found: {scores_file}")
        return None
    
    ranked_genes = load_ranked_genes(scores_file)
    print(f"Ranked {len(ranked_genes)} genes by effect score")
    
    # Test pathways
    results = []
    
    print(f"\n{'Pathway':<25} {'N':<6} {'Median Rank':<12} {'P-value':<12}")
    print("-" * 60)
    
    for pathway_name, pathway_genes in PATHWAY_SETS.items():
        stat, p_val, effect, med_rank, n_found = gsea_mannwhitney(
            ranked_genes, pathway_genes
        )
        
        if p_val is not None:
            sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
            print(f"{pathway_name:<25} {n_found:<6} {med_rank:<12.1f} {p_val:<12.4f} {sig}")
            
            results.append({
                'Pathway': pathway_name,
                'N_Found': n_found,
                'MedianRank': med_rank,
                'P_value': p_val,
                'RankBiserial': effect
            })
    
    results_df = pd.DataFrame(results)
    
    # FDR correction
    from statsmodels.stats.multitest import multipletests
    results_df['P_FDR'] = multipletests(results_df['P_value'], method='fdr_bh')[1]
    
    # Save
    output_file = TABLES_DIR / "gsea_results.csv"
    results_df.to_csv(output_file, index=False)
    print(f"\n✓ Saved to {output_file}")
    
    # Plot
    create_visualization(results_df)
    
    return results_df


def create_visualization(results_df):
    """Create lollipop plot of pathway enrichment."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    df = results_df.sort_values('P_value')
    y = range(len(df))
    x = -np.log10(df['P_value'])
    
    colors = ['#2ecc71' if p < 0.05 else '#95a5a6' for p in df['P_value']]
    
    ax.hlines(y, 0, x, colors=colors, linewidth=2)
    ax.scatter(x, y, c=colors, s=100, zorder=5)
    ax.axvline(-np.log10(0.05), color='red', linestyle='--', label='P = 0.05')
    
    ax.set_yticks(y)
    ax.set_yticklabels(df['Pathway'])
    ax.set_xlabel('-log₁₀(P-value)', fontsize=12)
    ax.set_title('Ranked GSEA: Pathway Enrichment in SCZ Risk Genes', fontsize=13)
    ax.legend()
    
    plt.tight_layout()
    output_file = FIGURES_DIR / "gsea_enrichment.png"
    plt.savefig(output_file, dpi=150)
    plt.close()
    print(f"✓ Figure saved to {output_file}")


if __name__ == "__main__":
    main()
