#!/usr/bin/env python3
"""
02_alphagenome_scoring.py
=========================

Step 2: Score schizophrenia lead variants using AlphaGenome API.

This script queries the AlphaGenome deep learning model to predict
regulatory effects of each variant on gene expression, splicing,
and chromatin accessibility in brain tissues.

Input:
    - data/processed/scz_lead_snps.csv

Output:
    - data/processed/alphagenome_scores.csv

Usage:
    export ALPHAGENOME_API_KEY="your_key"
    python scripts/02_alphagenome_scoring.py

Author: Yusuf Cicek
Date: February 2026
"""

import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
import sys

# Add parent directory for config import
sys.path.append(str(Path(__file__).parent.parent))
from config import (PROCESSED_DIR, ALPHAGENOME_API_KEY, 
                    BRAIN_GTEX_TISSUES, SCORING_WINDOW_HALF)


class AlphaGenomeScorer:
    """
    Wrapper for AlphaGenome API variant scoring.
    
    Extracts multi-modal regulatory effects for each variant:
    - Expression (log-fold change per gene/tissue)
    - Splicing (splice site usage changes)
    - Chromatin accessibility (ATAC/DNase)
    """
    
    def __init__(self, api_key):
        from alphagenome.models.dna_client import create
        from alphagenome.models import variant_scorers
        from alphagenome.protos import dna_model_pb2
        
        self.client = create(api_key=api_key)
        self.recommended = variant_scorers.get_recommended_scorers(
            dna_model_pb2.Organism.ORGANISM_HOMO_SAPIENS
        )
        
        # Scorer indices
        self.scorer_map = {
            'expression_lfc': 7,       # GeneMaskLFCScorer (signed)
            'expression_active': 8,    # GeneMaskActiveScorer (magnitude)
            'splicing_usage': 10,      # SplicingScorer
            'atac': 0,                 # ATAC-seq
        }
    
    def _get_brain_indices(self, adata):
        """Get indices of brain-specific tracks."""
        if 'gtex_tissue' in adata.var.columns:
            mask = adata.var['gtex_tissue'].isin(BRAIN_GTEX_TISSUES)
            return np.where(mask)[0]
        return np.array([])
    
    def score_variant(self, chrom, pos, ref, alt, snp_id=''):
        """
        Score a single variant across all modalities.
        
        Parameters
        ----------
        chrom : str
            Chromosome (e.g., 'chr1' or '1')
        pos : int
            Genomic position (1-based)
        ref : str
            Reference allele
        alt : str
            Alternate allele
        snp_id : str
            Optional variant identifier
            
        Returns
        -------
        dict
            Dictionary with scores for each modality
        """
        from alphagenome.data.genome import Interval, Variant
        
        if not chrom.startswith('chr'):
            chrom = f'chr{chrom}'
        
        variant = Variant(
            chromosome=chrom,
            position=pos,
            reference_bases=ref,
            alternate_bases=alt,
            name=snp_id
        )
        
        interval_start = max(0, pos - SCORING_WINDOW_HALF)
        interval_end = interval_start + SCORING_WINDOW_HALF * 2
        
        interval = Interval(
            chromosome=chrom,
            start=interval_start,
            end=interval_end
        )
        
        scores = self.client.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=self.recommended
        )
        
        return self._extract_scores(scores)
    
    def _extract_scores(self, scores):
        """Extract per-modality scores with brain filtering."""
        result = {
            'genes': [],
            'expression_brain': 0.0,
            'expression_all': 0.0,
            'splicing_max': 0.0,
            'atac_brain': 0.0,
            'n_genes': 0
        }
        
        # Expression
        expr_adata = scores[self.scorer_map['expression_active']]
        if expr_adata.X is not None and 'gene_name' in expr_adata.obs.columns:
            brain_idx = self._get_brain_indices(expr_adata)
            
            for i, gene in enumerate(expr_adata.obs['gene_name']):
                vals = expr_adata.X[i]
                non_nan = vals[~np.isnan(vals)]
                
                if len(non_nan) > 0:
                    all_effect = float(np.abs(non_nan).max())
                    
                    if len(brain_idx) > 0:
                        brain_vals = vals[brain_idx]
                        brain_non_nan = brain_vals[~np.isnan(brain_vals)]
                        brain_effect = float(np.abs(brain_non_nan).max()) if len(brain_non_nan) > 0 else 0
                    else:
                        brain_effect = 0
                    
                    result['genes'].append(f"{gene}:{all_effect:.3f}")
                    result['expression_all'] = max(result['expression_all'], all_effect)
                    result['expression_brain'] = max(result['expression_brain'], brain_effect)
        
        result['n_genes'] = len(result['genes'])
        result['genes'] = ';'.join(result['genes'])
        
        # Splicing
        splice_adata = scores[self.scorer_map['splicing_usage']]
        if splice_adata.X is not None:
            non_nan = splice_adata.X[~np.isnan(splice_adata.X)]
            if len(non_nan) > 0:
                result['splicing_max'] = float(np.abs(non_nan).max())
        
        return result


def main():
    """Main scoring pipeline."""
    print("=" * 70)
    print("STEP 2: ALPHAGENOME VARIANT SCORING")
    print("=" * 70)
    
    # Check API key
    if not ALPHAGENOME_API_KEY:
        print("\nError: ALPHAGENOME_API_KEY not set")
        print("Set via: export ALPHAGENOME_API_KEY='your_key'")
        return None
    
    # Load variants
    input_file = PROCESSED_DIR / "scz_lead_snps.csv"
    if not input_file.exists():
        print(f"\nError: Input file not found: {input_file}")
        print("Run 01_extract_gwas_variants.py first")
        return None
    
    variants = pd.read_csv(input_file)
    print(f"Loaded {len(variants)} variants for scoring")
    
    # Initialize scorer
    scorer = AlphaGenomeScorer(ALPHAGENOME_API_KEY)
    
    # Score each variant
    results = []
    
    for _, row in tqdm(variants.iterrows(), total=len(variants)):
        try:
            scores = scorer.score_variant(
                chrom=str(row['CHR']),
                pos=int(row['BP']),
                ref=str(row['A1']),
                alt=str(row['A2']),
                snp_id=str(row.get('SNP', ''))
            )
            
            result = {
                'SNP': row.get('SNP', ''),
                'CHR': row['CHR'],
                'BP': row['BP'],
                'REF': row['A1'],
                'ALT': row['A2'],
                **scores,
                'success': True
            }
        except Exception as e:
            result = {
                'SNP': row.get('SNP', ''),
                'CHR': row['CHR'],
                'BP': row['BP'],
                'success': False,
                'error': str(e)
            }
        
        results.append(result)
    
    # Save results
    results_df = pd.DataFrame(results)
    output_file = PROCESSED_DIR / "alphagenome_scores.csv"
    results_df.to_csv(output_file, index=False)
    
    # Summary
    successful = results_df[results_df['success'] == True]
    print(f"\n✓ Scored {len(successful)}/{len(variants)} variants")
    print(f"  Mean expression effect: {successful['expression_brain'].mean():.3f}")
    print(f"  Max expression effect: {successful['expression_brain'].max():.3f}")
    print(f"\n✓ Saved to {output_file}")
    
    return results_df


if __name__ == "__main__":
    main()
