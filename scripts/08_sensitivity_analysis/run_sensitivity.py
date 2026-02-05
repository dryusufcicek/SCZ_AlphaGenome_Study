
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import sys
from tqdm import tqdm

sys.path.append(str(Path(__file__).parent.parent.parent))
from config import PROCESSED_DIR, RESULTS_DIR

# Import GSEA logic from existing script
import importlib.util
spec = importlib.util.spec_from_file_location("gsea_module", Path(__file__).parent.parent.parent / "34_revision2_gsea_pathways_universe.py")
gsea_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gsea_module)

def run_sensitivity_analysis():
    print("=" * 60)
    print("SENSITIVITY ANALYSIS: Robustness to Outlier Removal")
    print("=" * 60)
    
    # 1. Load Z-Scores
    z_score_file = PROCESSED_DIR / "gene_z_scores.csv"
    if not z_score_file.exists():
        return
        
    df = pd.read_csv(z_score_file)
    
    # 2. Identify and Drop Top 1%
    threshold = df['z_score'].quantile(0.99)
    top_genes = df[df['z_score'] > threshold]['gene'].tolist()
    
    # Remove outliers from the HIT list
    df_hits_filtered = df[df['z_score'] <= threshold].copy()
    
    # 3. Construct the Full Universe (Hits + Zeros)
    # The original analysis used ~47,800 genes.
    # We will assume any gene in the GMT but NOT in our hits is a "Zero".
    
    # Load GMT
    gmt_file = PROCESSED_DIR / "GO_Biological_Process_2023.gmt"
    # Ensure it exists
    if not gmt_file.exists():
         gsea_module.download_gmt_if_needed()
         
    pathways = gsea_module.load_gmt(gmt_file)
    
    # Collect all unique genes in the GMT universe (Background)
    all_pathway_genes = set()
    for gs in pathways.values():
        all_pathway_genes.update(gs)
    
    # Map Hit Scores
    gene_scores = dict(zip(df_hits_filtered['gene'], df_hits_filtered['z_score']))
    
    # Build Full Rank List
    # Any gene in the GMT that isn't a hit gets 0.0
    # Any gene in our filtered hits keeps its score
    full_universe_scores = []
    
    # We consider the universe to be (GMT Genes UNION Hit Genes)
    full_universe_set = all_pathway_genes.union(set(gene_scores.keys()))
    
    for gene in full_universe_set:
        full_universe_scores.append({
            'gene': gene,
            'z_score': gene_scores.get(gene, 0.0) # 0.0 if not in hits
        })
        
    # Convert to DataFrame for Ranking
    rank_df = pd.DataFrame(full_universe_scores)
    
    # Assign Ranks (Higher Z = Better Rank = Lower Number)
    # MannWhitneyU in the GSEA script tested 'less', meaning lower rank is better.
    # So we sort descending by score.
    rank_df = rank_df.sort_values('z_score', ascending=False).reset_index(drop=True)
    rank_df['rank'] = rank_df.index + 1
    
    print(f"Full Universe Size: {len(rank_df)}")
    print(f"Non-Zero Hits Remaining: {len(df_hits_filtered)}")
    
    # 4. Re-Run GSEA
    print("\nRunning GSEA on Filtered Data...")
    results = []
    
    for name, genes in tqdm(pathways.items()):
        # Call the imported function
        p, r, med, n = gsea_module.gsea_mannwhitney(rank_df, genes)
        
        if p is not None:
             results.append((name, p))
             
    results.sort(key=lambda x: x[1])
    
    # 5. Report
    print("\n" + "=" * 60)
    print("RESULTS: Top Pathways (Top 1% Removed)")
    print("=" * 60)
    for i, (name, p) in enumerate(results[:20]):
        print(f"{i+1}. {name[:50]} (P = {p:.2e})")

if __name__ == "__main__":
    run_sensitivity_analysis()
