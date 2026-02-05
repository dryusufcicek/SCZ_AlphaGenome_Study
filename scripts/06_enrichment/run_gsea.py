#!/usr/bin/env python3
"""
34_revision2_gsea_pathways_universe.py
======================================

Revision 2 (Fix): Ranked GSEA with WHOLE GENOME UNIVERSE & UNBIASED GO DISCOVERY.

Addressing Reviewer Critiques:
1. "Universe Bias": Uses 47k gene zero-inflated background.
2. "Cherry Picking": Tests ALL ~5,000 GO Biological Processes (Unbiased).

Methodology:
1. Load "Universe" (~20k - 47k genes).
2. Load AlphaGenome scores (Detected=Score, Others=0.0).
3. Download/Load GO_Biological_Process_2023.gmt.
4. Perform Mann-Whitney U on ALL gene sets.
5. Apply FDR correction across the full discovery set.

"""

import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import sys
import requests
from tqdm import tqdm

sys.path.append(str(Path(__file__).parent.parent.parent))
from config import PROCESSED_DIR, RESULTS_DIR

# URL for Gene Ontology (Enrichr / MaayanLab)
GO_URL = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2023"
GMT_FILE = PROCESSED_DIR / "GO_Biological_Process_2023.gmt"

def download_gmt_if_needed():
    """Download GO database if not present."""
    if GMT_FILE.exists():
        print(f"Using existing GMT: {GMT_FILE}")
        return

    print(f"Downloading GO Biological Process from {GO_URL}...")
    try:
        r = requests.get(GO_URL)
        r.raise_for_status()
        with open(GMT_FILE, 'w') as f:
            f.write(r.text)
        print("✅ Download complete.")
    except Exception as e:
        print(f"❌ Failed to download GMT: {e}")
        sys.exit(1)

def load_gmt(gmt_path):
    """Parse GMT file into dict."""
    pathways = {}
    print(f"Parsing {gmt_path}...")
    with open(gmt_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3: continue
            name = parts[0]
            # Skip description (parts[1])
            genes = [g.split(',')[0] for g in parts[2:] if g] # Handle occasional gene,1.0 format
            
            # Filter for reasonable size (avoid overly broad or specific terms)
            if 15 <= len(genes) <= 500:
                pathways[name] = genes
                
    print(f"Loaded {len(pathways)} pathways (filtered 15-500 genes).")
    return pathways

def load_universe_genes():
    """Load list of all considered genes (approx 20k)."""
    # Try multiple sources including JSI and PGC3 lists
    sources = [
        PROCESSED_DIR / "gene_jsi.csv",
        PROCESSED_DIR / "pgc3_all_genes.csv",
        PROCESSED_DIR / "gene_developmental_peaks.csv"
    ]
    
    universe = set()
    for s in sources:
        if s.exists():
            df = pd.read_csv(s)
            col = next((c for c in df.columns if 'gene' in c.lower() or 'symbol' in c.lower()), None)
            if col:
                universe.update(df[col].dropna().unique())
                
    print(f"Loaded Universe: {len(universe)} genes")
    return list(universe)

def load_ranked_genes_with_universe(scores_file, universe_list):
    """Load scores and map to UNIVERSE with 0-inflation."""
    if not scores_file.exists():
        print(f"Error: Scores file not found: {scores_file}")
        return None
    
    # Parse Scores
    gene_scores = {}
    try:
        scores = pd.read_csv(scores_file)
    except:
        return None

    # Handle both column formats
    gene_col = 'all_genes' if 'all_genes' in scores.columns else 'genes'
    
    for _, row in tqdm(scores.iterrows(), total=len(scores), desc="Parsing Scores"):
        val = row.get(gene_col)
        if pd.notna(val) and val:
            for entry in str(val).split(';'):
                if ':' in entry:
                    parts = entry.split(':')
                    try:
                        g = parts[0]
                        s = float(parts[1])
                        # Keep max score
                        gene_scores[g] = max(gene_scores.get(g, 0), s)
                    except:
                        pass

    # Create Full Ranking
    full_list = []
    
    # Ensure universe includes all detected genes too
    final_universe = set(universe_list) | set(gene_scores.keys())
    print(f"Final Universe Size for Ranking: {len(final_universe)}")
    
    for gene in final_universe:
        score = gene_scores.get(gene, 0.0)
        full_list.append({'gene': gene, 'score': score})
        
    df = pd.DataFrame(full_list).sort_values('score', ascending=False).reset_index(drop=True)
    df['rank'] = range(1, len(df) + 1)
    
    return df

def gsea_mannwhitney(ranked_genes, pathway_genes):
    """Mann-Whitney U Test."""
    pathway_set = set(pathway_genes)
    
    # Intersection with our universe
    valid_genes = [g for g in pathway_genes if g in ranked_genes['gene'].values]
    
    if len(valid_genes) < 5: # SKIP tiny overlap
        return None, None, 0, len(valid_genes)
        
    pathway_ranks = ranked_genes[ranked_genes['gene'].isin(pathway_set)]['rank'].values
    non_pathway_ranks = ranked_genes[~ranked_genes['gene'].isin(pathway_set)]['rank'].values
    
    # Test: Pathway ranks are LOWER (better) than background?
    stat, p_val = stats.mannwhitneyu(pathway_ranks, non_pathway_ranks, alternative='less')
    
    # Effect Size (Rank Biserial)
    n1, n2 = len(pathway_ranks), len(non_pathway_ranks)
    r = 1 - (2 * stat) / (n1 * n2)
    
    return p_val, r, np.median(pathway_ranks), len(valid_genes)

def main():
    print("=" * 70)
    print("STEP 3: UNBIASED GSEA (WHOLE GENOME + FULL GO DATABASE)")
    print("=" * 70)
    
    # 1. Prepare Data
    download_gmt_if_needed()
    universe = load_universe_genes()
    
    scores_file = PROCESSED_DIR / "alphgenome_variant_scores_proxy.csv"
    if not scores_file.exists():
        print("⚠️ Proxy scores not ready. Utilizing partial file if available.")
        scores_file = PROCESSED_DIR / "alphgenome_variant_scores.csv"

    ranked_genes = load_ranked_genes_with_universe(scores_file, universe)
    if ranked_genes is None: return

    pathways = load_gmt(GMT_FILE)
    
    # 2. Run Unbiased Scan
    print(f"\nScanning {len(pathways)} pathways against {len(ranked_genes)} genes...")
    results = []
    
    # Progress bar
    for name, genes in tqdm(pathways.items(), desc="Running GSEA"):
        p, r, med, n = gsea_mannwhitney(ranked_genes, genes)
        if p is not None:
             results.append({
                 'Pathway': name,
                 'N_Overlap': n,
                 'P_Value': p,
                 'Effect_Size': r,
                 'Median_Rank': med
             })
             
    # 3. Process Results
    res_df = pd.DataFrame(results)
    
    if not res_df.empty:
        # FDR Correction
        from statsmodels.stats.multitest import multipletests
        res_df['FDR'] = multipletests(res_df['P_Value'], method='fdr_bh')[1]
        
        # Sort by P-value
        res_df = res_df.sort_values('P_Value')
        
        # Display Top 50
        print("\n" + "="*80)
        print("TOP 20 ENRICHED PATHWAYS (UNBIASED DISCOVERY)")
        print("="*80)
        print(f"{'Pathway':<50} {'N':<5} {'P-Value':<10} {'FDR':<10} {'Effect':<6}")
        print("-" * 80)
        
        for _, row in res_df.head(20).iterrows():
            clean_name = row['Pathway'].split('(GO')[0][:48]
            print(f"{clean_name:<50} {row['N_Overlap']:<5} {row['P_Value']:<10.2e} {row['FDR']:<10.2e} {row['Effect_Size']:<6.2f}")
            
        print("-" * 80)
        
        n_sig = len(res_df[res_df['FDR'] < 0.05])
        print(f"\nTotal Significant Pathways (FDR < 0.05): {n_sig}")
        
        # Save
        out = RESULTS_DIR / "gsea_unbiased_results.csv"
        res_df.to_csv(out, index=False)
        print(f"\nSaved full results to {out}")

if __name__ == "__main__":
    main()
