
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import sys
from tqdm import tqdm

sys.path.append(str(Path(__file__).parent.parent))
from config import PROCESSED_DIR, RESULTS_DIR, GENE_DIR

# Import GSEA logic from the new location
import importlib.util
gsea_script = Path(__file__).parent.parent / "06_enrichment" / "run_gsea.py"
spec = importlib.util.spec_from_file_location("gsea_module", gsea_script)
gsea_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gsea_module)

def run_sensitivity_analysis(exclude_top_1pct=True, exclude_mhc=False):
    title = "SENSITIVITY ANALYSIS"
    if exclude_top_1pct: title += ": Outlier Removal (Top 1%)"
    if exclude_mhc: title += " + MHC Exclusion (Chr6:25-34Mb)"
    
    print("\n" + "=" * 60)
    print(title)
    print("=" * 60)
    
    # 1. Load Z-Scores
    z_score_file = GENE_DIR / "gene_z_scores.csv"
    if not z_score_file.exists():
        print(f"Error: {z_score_file} not found.")
        return
        
    df = pd.read_csv(z_score_file)

    
    # 2. Filter Genes
    filtered_df = df.copy()
    
    if exclude_top_1pct:
        threshold = df['z_score'].quantile(0.99)
        filtered_df = filtered_df[filtered_df['z_score'] <= threshold]
        print(f"Dropped Top 1% genes (Z > {threshold:.2f})")
        
    if exclude_mhc:
        # We need chromosome info. Re-map if needed or use gene names
        # Standard MHC genes often start with HLA-, HIST*, etc.
        # But for rigor, we should use genomic coordinates if possible.
        # Since we don't have CHR in gene_z_scores.csv, we rely on a mapping file.
        mapping_file = PROCESSED_DIR / "variant_gene_mapping.csv"
        if mapping_file.exists():
            mapping = pd.read_csv(mapping_file)
            mhc_genes = mapping[mapping['CHR'].isin(['6', 'chr6']) & (mapping['BP'] >= 25000000) & (mapping['BP'] <= 34000000)]['gene'].unique()
            filtered_df = filtered_df[~filtered_df['gene'].isin(mhc_genes)]
            print(f"Excluded {len(mhc_genes)} genes from MHC region.")
        else:
            # Fallback simple string match for Histones in MHC
            mhc_pattern = "^HLA-|^H1|^H2|^H3|^H4"
            mhc_genes = filtered_df[filtered_df['gene'].str.contains(mhc_pattern)]['gene'].unique()
            filtered_df = filtered_df[~filtered_df['gene'].isin(mhc_genes)]
            print(f"Excluded {len(mhc_genes)} genes (Histone/HLA proxy) due to missing mapping file.")

    # 3. Construct the Full Universe
    gmt_file = PROCESSED_DIR / "GO_Biological_Process_2023.gmt"
    if not gmt_file.exists():
         gsea_module.download_gmt_if_needed()
         
    pathways = gsea_module.load_gmt(gmt_file)
    
    all_pathway_genes = set()
    for gs in pathways.values():
        all_pathway_genes.update(gs)
    
    gene_scores = dict(zip(filtered_df['gene'], filtered_df['z_score']))
    full_universe_set = all_pathway_genes.union(set(gene_scores.keys()))
    
    full_universe_scores = []
    for gene in full_universe_set:
        full_universe_scores.append({
            'gene': gene,
            'z_score': gene_scores.get(gene, 0.0)
        })
        
    rank_df = pd.DataFrame(full_universe_scores)
    rank_df = rank_df.sort_values('z_score', ascending=False).reset_index(drop=True)
    rank_df['rank'] = rank_df.index + 1
    
    # 4. Re-Run GSEA
    results = []
    for name, genes in tqdm(pathways.items(), desc="Running GSEA"):
        p, r, med, n = gsea_module.gsea_mannwhitney(rank_df, genes)
        if p is not None:
             results.append((name, p))
             
    results.sort(key=lambda x: x[1])
    
    # 5. Save and Report
    suffix = "_no_top1pct" if exclude_top_1pct else ""
    if exclude_mhc: suffix += "_no_mhc"
    
    out_file = RESULTS_DIR / f"sensitivity_results{suffix}.csv"
    pd.DataFrame(results, columns=['Pathway', 'P-value']).to_csv(out_file, index=False)
    
    print(f"\nRESULTS: Top Pathways {suffix.replace('_', ' ')}")
    for i, (name, p) in enumerate(results[:15]):
        print(f"{i+1}. {name[:50]} (P = {p:.2e})")

def main():
    # Run the core PI-requested sensitivity tests
    run_sensitivity_analysis(exclude_top_1pct=True, exclude_mhc=False)
    run_sensitivity_analysis(exclude_top_1pct=False, exclude_mhc=True)
    run_sensitivity_analysis(exclude_top_1pct=True, exclude_mhc=True)

if __name__ == "__main__":
    main()

