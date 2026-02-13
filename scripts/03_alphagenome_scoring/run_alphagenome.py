"""
Step 3: Score SCZ Fine-Mapped Variants using AlphaGenome API

Input: Official PGC3 FINEMAP 95% credible sets (20,591 variants)
Output: Variant-level regulatory effect scores with posterior probabilities (PP)
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
    Returns DataFrame with results.
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
    output_file = PROCESSED_DIR / "alphagenome_variant_scores.csv"
    seen_snps = set()
    if output_file.exists():
        try:
            existing = pd.read_csv(output_file)
            seen_snps = set(existing['SNP'].astype(str))
            print(f"Resuming: Found {len(seen_snps)} already scored variants.")
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
            
            # Extract Max Gene Effect & Modalities
            target_gene = None
            target_gene_score = 0.0
            all_genes_str = ""
            track_scores = {t: 0.0 for t in ["DNase", "H3K27ac", "H3K4me1", "H3K4me3", "CAGE", "CTCF", "RNA"]}
            
            if gene_scorer_idx is not None and len(scores) > gene_scorer_idx:
                gene_adata = scores[gene_scorer_idx]
                if gene_adata.X is not None and 'gene_name' in gene_adata.obs.columns:
                    gene_effects = []
                    track_names = ["DNase", "H3K27ac", "H3K4me1", "H3K4me3", "CAGE", "CTCF", "RNA"]
                    for j, gname in enumerate(gene_adata.obs['gene_name']):
                        vals = gene_adata.X[j]  # These are the 7 tracks
                        
                        # Calculate robust max for the gene
                        valid_vals = vals[~np.isnan(vals)]
                        if len(valid_vals) > 0:
                            gene_max_effect = float(np.abs(valid_vals).max())
                            
                            # Keep track for all genes string (max-based)
                            if gene_max_effect > 0.1:
                                gene_effects.append(f"{gname}:{gene_max_effect:.3f}")
                            
                            # Update target gene if this gene has the global max
                            if gene_max_effect > target_gene_score:
                                target_gene_score = gene_max_effect
                                target_gene = gname
                                # Save individual tracks for the champion gene
                                for t_idx, t_name in enumerate(track_names):
                                    if t_idx < len(vals):
                                        track_scores[t_name] = float(vals[t_idx])
                    
                    all_genes_str = ";".join(gene_effects)

            # Record Result
            res = {
                'SNP': snp_id,
                'target_gene': target_gene,
                'target_gene_score': target_gene_score,
                'all_genes': all_genes_str,
                'success': True
            }
            # Add track columns
            for t_name, t_val in track_scores.items():
                res[f'score_{t_name}'] = t_val
            
            temp_results.append(res)
            
        except Exception as e:
            # Log failure but continue
            print(f"Error scoring {snp_id}: {e}")
            temp_results.append({
                'SNP': snp_id,
                'success': False, 
                'error': str(e)
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
    print("Step 3: AlphaGenome Scoring (Official PGC3 Credible Sets)")
    print("="*60)
    
    # 1. Load Official PGC3 Credible Sets (Target list)
    official_file = PROCESSED_DIR / "credible_sets" / "scz_credible_sets_official.csv"
    if not official_file.exists():
        print(f"Input file not found: {official_file}")
        return
    official_df = pd.read_csv(official_file)
    print(f"Official Target: {len(official_df)} variants")
    
    # 2. Score ALL variants (High-Res 7-Modalities)
    print("Starting FRESH run for ALL variants (7-Modalities)...")
    
    api_key = ALPHAGENOME_API_KEY or os.environ.get('ALPHAGENOME_API_KEY')
    if not api_key:
        print("❌ API Key missing.")
        return

    # Run Scoring
    scored_df = score_variants_with_alphagenome(official_df, api_key)
    
    # 3. Merge with Metadata (PP) and Save
    print("\nMerging metadata and saving...")
    
    # Drop columns that will be added from official_df to avoid duplication
    cols_to_drop = ['CHR', 'BP', 'A1', 'A2', 'Locus_ID', 'PP'] 
    scored_df = scored_df.drop(columns=[c for c in cols_to_drop if c in scored_df.columns], errors='ignore')
    
    # Merge
    final_df = official_df.merge(scored_df, on='SNP', how='inner')
    
    out_file = PROCESSED_DIR / "alphagenome_variant_scores.csv"
    final_df.to_csv(out_file, index=False)
    print(f"\n✅ Saved {len(final_df)} variants to {out_file}")
    print("Full 7-modality dataset generated.")

if __name__ == "__main__":
    main()
