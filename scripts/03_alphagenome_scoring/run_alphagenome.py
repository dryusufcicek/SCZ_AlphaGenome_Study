"""
Step 3 (REVISED): Score SCZ "Proxy Credible Set" variants using AlphaGenome API

Updates:
- Consumes `scz_credible_sets_proxy.csv` (10k+ variants) instead of Lead SNPs
- Preserves `Proxy_PP` and `Locus_ID` for downstream weighted aggregation
- Handles larger batch size with progress saving
"""
import pandas as pd
import numpy as np
from pathlib import Path
import sys
import os
import time
from tqdm import tqdm

sys.path.append(str(Path(__file__).parent.parent.parent))
from config import PROCESSED_DIR, ALPHAGENOME_API_KEY

def score_variants_with_alphagenome(variants_df, api_key):
    """
    Score variants using AlphaGenome API.
    """
    from alphagenome.models.dna_client import create
    from alphagenome.data.genome import Interval, Variant
    from alphagenome.models import variant_scorers
    from alphagenome.protos import dna_model_pb2

    print(f"Initializing AlphaGenome client...")
    try:
        client = create(api_key=api_key)
    except Exception as e:
        print(f"Failed to create client: {e}")
        return pd.DataFrame()

    # Get recommended scorers
    recommended = variant_scorers.get_recommended_scorers(
        dna_model_pb2.Organism.ORGANISM_HOMO_SAPIENS
    )
    
    # Find GeneMaskActiveScorer
    gene_scorer_idx = next((i for i, s in enumerate(recommended) if type(s).__name__ == 'GeneMaskActiveScorer'), None)
    
    results = []
    
    # Check for existing partial results to resume
    output_file = PROCESSED_DIR / "alphgenome_variant_scores_proxy.csv"
    seen_snps = set()
    if output_file.exists():
        try:
            existing = pd.read_csv(output_file)
            seen_snps = set(existing['SNP'].astype(str))
            print(f"Resuming: Found {len(seen_snps)} already scored variants.")
            # We will append new results
            # Note: valid restart logic is complex, simpler to just skip in loop
        except:
            pass

    print(f"\nScoring {len(variants_df)} variants...")
    
    # Batch save every 100 variants
    BATCH_SIZE = 100
    temp_results = []
    
    for idx, row in tqdm(variants_df.iterrows(), total=len(variants_df), desc="Scoring"):
        snp_id = str(row['SNP'])
        
        # Skip if already done
        if snp_id in seen_snps:
            continue
            
        try:
            # Prepare Inputs
            chrom = str(row['CHR'])
            if not chrom.startswith('chr'): chrom = f'chr{chrom}'
            pos = int(row['BP'])
            ref = str(row['A1'])
            alt = str(row['A2'])
            
            variant = Variant(chromosome=chrom, position=pos, reference_bases=ref, alternate_bases=alt, name=snp_id)
            
            # Window 512kb
            interval_start = max(0, pos - 262144)
            interval_end = interval_start + 524288
            interval = Interval(chromosome=chrom, start=interval_start, end=interval_end)
            
            # API Call
            scores = client.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=recommended
            )
            
            # Extract Max Gene Effect
            target_gene = None
            target_gene_score = 0.0
            all_genes_str = ""
            
            if gene_scorer_idx is not None and len(scores) > gene_scorer_idx:
                gene_adata = scores[gene_scorer_idx]
                if gene_adata.X is not None and 'gene_name' in gene_adata.obs.columns:
                    gene_effects = []
                    for j, gname in enumerate(gene_adata.obs['gene_name']):
                        vals = gene_adata.X[j]
                        # robust max
                        valid_vals = vals[~np.isnan(vals)]
                        if len(valid_vals) > 0:
                            effect = float(np.abs(valid_vals).max())
                            # Threshold optimization to save space
                            if effect > 0.1: 
                                gene_effects.append(f"{gname}:{effect:.3f}")
                                if effect > target_gene_score:
                                    target_gene_score = effect
                                    target_gene = gname
                    all_genes_str = ";".join(gene_effects)

            # Record Result
            res = {
                'SNP': snp_id,
                'CHR': row['CHR'],
                'BP': pos,
                'A1': ref, 
                'A2': alt,
                'Locus_ID': row.get('Locus_ID', -1),
                'Proxy_PP': row.get('Proxy_PP', 0), # Critical for weighting
                'target_gene': target_gene,
                'target_gene_score': target_gene_score,
                'all_genes': all_genes_str,
                'success': True
            }
            temp_results.append(res)
            
        except Exception as e:
            # Log failure but continue
            temp_results.append({
                'SNP': snp_id, 'CHR': row['CHR'], 'BP': pos,
                'Locus_ID': row.get('Locus_ID', -1),
                'Proxy_PP': row.get('Proxy_PP', 0),
                'success': False, 'error': str(e)
            })
            
        # Incremental Save
        if len(temp_results) >= BATCH_SIZE:
            batch_df = pd.DataFrame(temp_results)
            # Append to CSV
            hdr = not output_file.exists()
            batch_df.to_csv(output_file, mode='a', header=hdr, index=False)
            temp_results = []
            
    # Final Save
    if temp_results:
        batch_df = pd.DataFrame(temp_results)
        hdr = not output_file.exists()
        batch_df.to_csv(output_file, mode='a', header=hdr, index=False)

    print(f"\nCompleted. Results saved to {output_file}")
    return pd.read_csv(output_file)

def main():
    print("="*60)
    print("Step 3: AlphaGenome Scoring (Proxy Credible Sets)")
    print("="*60)
    
    # Load Input
    input_file = PROCESSED_DIR / "scz_credible_sets_proxy.csv"
    if not input_file.exists():
        print(f"Input file not found: {input_file}")
        print("Please run Step 1 (01_extract_gwas_variants.py) first.")
        return

    variants = pd.read_csv(input_file)
    print(f"Loaded {len(variants)} variants from Proxy Credible Sets")
    
    # API Check
    api_key = ALPHAGENOME_API_KEY or os.environ.get('ALPHAGENOME_API_KEY')
    if not api_key:
        print("❌ API Key missing.")
        return

    # Run
    score_variants_with_alphagenome(variants, api_key)

if __name__ == "__main__":
    main()
