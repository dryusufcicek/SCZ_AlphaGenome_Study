"""
Revision 3b: Signed Score Extraction for Directionality Analysis

CRITICAL FIX: The original multimodal_scoring.py used np.abs() to store effect
magnitudes, discarding sign information. This script re-extracts SIGNED scores
to properly characterize:
- Positive scores: Risk allele → INCREASED expression (Gain-of-Function)
- Negative scores: Risk allele → DECREASED expression (Loss-of-Function)

This is essential for determining drug target directionality (agonist vs antagonist).
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import os
from tqdm import tqdm
from scipy import stats
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).parent))
from config import PROCESSED_DIR, RESULTS_DIR, ALPHAGENOME_API_KEY

# Brain tissue keywords
BRAIN_GTEX_TISSUES = [
    'Brain_Amygdala', 'Brain_Anterior_cingulate_cortex_BA24',
    'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere',
    'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9',
    'Brain_Hippocampus', 'Brain_Hypothalamus',
    'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Putamen_basal_ganglia',
    'Brain_Spinal_cord_cervical_c-1', 'Brain_Substantia_nigra'
]


class SignedScorer:
    """Extract SIGNED scores preserving effect direction."""
    
    def __init__(self, api_key):
        from alphagenome.models.dna_client import create
        from alphagenome.models import variant_scorers
        from alphagenome.protos import dna_model_pb2
        
        self.client = create(api_key=api_key)
        self.recommended = variant_scorers.get_recommended_scorers(
            dna_model_pb2.Organism.ORGANISM_HOMO_SAPIENS
        )
        
        # CRITICAL: Use LFC scorer for signed scores (not Active scorer)
        # Index 7: GeneMaskLFCScorer - provides log-fold-change (SIGNED!)
        # Index 8: GeneMaskActiveScorer - provides magnitude only (always positive)
        self.scorer_map = {
            'expression_lfc': 7,      # GeneMaskLFCScorer (log fold change) - SIGNED!
            'expression_active': 8,   # GeneMaskActiveScorer (magnitude only)
        }
    
    def _get_brain_track_indices(self, adata):
        """Get indices of brain-specific tracks."""
        if 'gtex_tissue' in adata.var.columns:
            brain_mask = adata.var['gtex_tissue'].isin(BRAIN_GTEX_TISSUES)
            return np.where(brain_mask)[0]
        return np.array([])
    
    def score_variant_signed(self, chrom, pos, ref, alt, snp_id=''):
        """Score variant preserving sign information."""
        from alphagenome.data.genome import Interval, Variant
        
        if not chrom.startswith('chr'):
            chrom = f'chr{chrom}'
        
        variant = Variant(
            chromosome=chrom, position=pos,
            reference_bases=ref, alternate_bases=alt, name=snp_id
        )
        
        window_half = 262144
        interval_start = max(0, pos - window_half)
        interval_end = interval_start + 524288
        
        interval = Interval(
            chromosome=chrom, start=interval_start, end=interval_end
        )
        
        scores = self.client.score_variant(
            interval=interval, variant=variant,
            variant_scorers=self.recommended
        )
        
        return self._extract_signed_scores(scores)
    
    def _extract_signed_scores(self, scores):
        """Extract SIGNED per-gene scores using LFC scorer."""
        result = {
            'genes': [],
            'signed_scores': [],
            'n_positive': 0,  # UP-regulated
            'n_negative': 0,  # DOWN-regulated
            'mean_signed': 0.0,
            'max_positive': 0.0,
            'min_negative': 0.0
        }
        
        # CRITICAL: Use LFC scorer (index 7) for signed scores
        # This returns log-fold-change with negative = decreased expression
        expr_adata = scores[self.scorer_map['expression_lfc']]
        
        if expr_adata.X is not None and 'gene_name' in expr_adata.obs.columns:
            brain_idx = self._get_brain_track_indices(expr_adata)
            
            for i, gene in enumerate(expr_adata.obs['gene_name']):
                vals = expr_adata.X[i]
                
                # Get brain-specific values
                if len(brain_idx) > 0:
                    brain_vals = vals[brain_idx]
                    brain_non_nan = brain_vals[~np.isnan(brain_vals)]
                else:
                    brain_non_nan = vals[~np.isnan(vals)]
                
                if len(brain_non_nan) > 0:
                    # KEEP SIGN! Use mean to preserve direction
                    signed_effect = float(np.mean(brain_non_nan))
                    
                    # Also track the extreme values (signed)
                    max_val = float(np.max(brain_non_nan))
                    min_val = float(np.min(brain_non_nan))
                    
                    # Use the extreme with larger absolute value
                    if abs(max_val) > abs(min_val):
                        dominant_effect = max_val
                    else:
                        dominant_effect = min_val
                    
                    result['genes'].append(gene)
                    result['signed_scores'].append(dominant_effect)
        
        # Summarize
        if result['signed_scores']:
            scores_arr = np.array(result['signed_scores'])
            result['n_positive'] = (scores_arr > 0).sum()
            result['n_negative'] = (scores_arr < 0).sum()
            result['mean_signed'] = float(np.mean(scores_arr))
            result['max_positive'] = float(np.max(scores_arr)) if result['n_positive'] > 0 else 0
            result['min_negative'] = float(np.min(scores_arr)) if result['n_negative'] > 0 else 0
        
        return result


def score_variants_signed(variants_df, api_key, max_variants=50):
    """Score variants with signed scores (limited for API efficiency)."""
    scorer = SignedScorer(api_key)
    all_genes = []
    all_scores = []
    variant_summary = []
    
    n_to_score = min(len(variants_df), max_variants)
    print(f"Scoring {n_to_score} variants with SIGNED scores...")
    
    for idx, row in tqdm(variants_df.head(n_to_score).iterrows(), total=n_to_score):
        try:
            chrom = str(row['CHR'])
            pos = int(row['BP'])
            ref = str(row.get('REF', row.get('A1', 'N')))
            alt = str(row.get('ALT', row.get('A2', 'N')))
            snp_id = str(row.get('SNP', f'{chrom}:{pos}'))
            
            signed_result = scorer.score_variant_signed(chrom, pos, ref, alt, snp_id)
            
            # Collect gene-level data
            for gene, score in zip(signed_result['genes'], signed_result['signed_scores']):
                all_genes.append(gene)
                all_scores.append(score)
            
            variant_summary.append({
                'SNP': snp_id,
                'n_genes': len(signed_result['genes']),
                'n_up': signed_result['n_positive'],
                'n_down': signed_result['n_negative'],
                'mean_signed': signed_result['mean_signed'],
                'success': True
            })
            
        except Exception as e:
            variant_summary.append({
                'SNP': str(row.get('SNP', 'unknown')),
                'n_genes': 0, 'n_up': 0, 'n_down': 0,
                'mean_signed': 0, 'success': False,
                'error': str(e)
            })
    
    return pd.DataFrame(variant_summary), all_genes, all_scores


def analyze_directionality(genes, scores):
    """Analyze LoF vs GoF patterns from signed scores."""
    print("\n" + "=" * 60)
    print("SIGNED DIRECTIONALITY ANALYSIS")
    print("=" * 60)
    
    scores_arr = np.array(scores)
    
    n_positive = (scores_arr > 0.1).sum()  # Meaningful UP effect
    n_negative = (scores_arr < -0.1).sum()  # Meaningful DOWN effect
    n_neutral = len(scores_arr) - n_positive - n_negative
    
    print(f"\nGene-level effect distribution (threshold ±0.1):")
    print(f"  UP (GoF):     {n_positive:>5} ({n_positive/len(scores_arr)*100:.1f}%)")
    print(f"  DOWN (LoF):   {n_negative:>5} ({n_negative/len(scores_arr)*100:.1f}%)")
    print(f"  NEUTRAL:      {n_neutral:>5} ({n_neutral/len(scores_arr)*100:.1f}%)")
    
    # Statistical test
    if n_positive + n_negative > 0:
        binom_p = stats.binomtest(n_negative, n_positive + n_negative, 0.5).pvalue
        print(f"\nBinomial test (deviation from 50/50): P = {binom_p:.2e}")
        
        if binom_p < 0.05:
            dominant = "DOWN (LoF)" if n_negative > n_positive else "UP (GoF)"
            print(f"✅ Significant directional bias toward {dominant}")
        else:
            print("→ No significant directional bias")
    
    # Mean effect direction
    print(f"\nOverall statistics:")
    print(f"  Mean signed effect: {np.mean(scores_arr):.4f}")
    print(f"  Std: {np.std(scores_arr):.4f}")
    
    # One-sample t-test: is mean significantly different from 0?
    t_stat, t_pval = stats.ttest_1samp(scores_arr, 0)
    print(f"  One-sample t-test (mean ≠ 0): t = {t_stat:.2f}, P = {t_pval:.2e}")
    
    return {
        'n_positive': n_positive,
        'n_negative': n_negative,
        'n_neutral': n_neutral,
        'mean_effect': np.mean(scores_arr),
        'binom_p': binom_p if n_positive + n_negative > 0 else None,
        't_test_p': t_pval
    }


def create_visualization(scores, output_path):
    """Create histogram of signed scores."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    scores_arr = np.array(scores)
    
    # Histogram with colored regions
    ax.hist(scores_arr, bins=50, color='#3498db', edgecolor='black', alpha=0.7)
    
    ax.axvline(x=0, color='red', linestyle='-', linewidth=2, label='Zero effect')
    ax.axvline(x=0.1, color='green', linestyle='--', linewidth=1.5, label='UP threshold')
    ax.axvline(x=-0.1, color='orange', linestyle='--', linewidth=1.5, label='DOWN threshold')
    ax.axvline(x=np.mean(scores_arr), color='purple', linestyle='-', linewidth=2, 
               label=f'Mean = {np.mean(scores_arr):.3f}')
    
    ax.set_xlabel('Signed Effect Score (+ = UP/GoF, - = DOWN/LoF)', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('Distribution of SIGNED AlphaGenome Effect Scores\n(Corrected Directionality Analysis)', fontsize=13)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"Visualization saved to {output_path}")


def main():
    print("=" * 70)
    print("REVISION 3b: Signed Score Extraction")
    print("=" * 70)
    
    # Check API key
    api_key = ALPHAGENOME_API_KEY or os.environ.get('ALPHAGENOME_API_KEY')
    if not api_key:
        print("ERROR: ALPHAGENOME_API_KEY not set")
        print("Please set the API key in config.py or environment variable")
        return
    
    # Load variants
    lead_snps_file = PROCESSED_DIR / "scz_lead_snps_robust.csv"
    if not lead_snps_file.exists():
        print(f"Lead SNPs file not found: {lead_snps_file}")
        return
    
    variants = pd.read_csv(lead_snps_file)
    print(f"Loaded {len(variants)} SCZ lead variants")
    
    # Score with signed extraction (limit to 50 for API efficiency)
    variant_results, all_genes, all_scores = score_variants_signed(
        variants, api_key, max_variants=50
    )
    
    if len(all_scores) == 0:
        print("No scores extracted. Check API connectivity.")
        return
    
    # Analyze directionality
    direction_stats = analyze_directionality(all_genes, all_scores)
    
    # Summary
    print("\n" + "=" * 70)
    print("CORRECTED DIRECTIONALITY FINDING")
    print("=" * 70)
    
    if direction_stats['n_negative'] > direction_stats['n_positive']:
        print("\n✅ SCZ risk alleles show LOSS-OF-FUNCTION bias")
        print("   This is consistent with SCHEMA rare variant findings")
    elif direction_stats['n_positive'] > direction_stats['n_negative']:
        print("\n✅ SCZ risk alleles show GAIN-OF-FUNCTION bias")
        print("   This suggests excitotoxicity / hyperactivity mechanism")
    else:
        print("\n→ Balanced directional effects (both LoF and GoF)")
        print("   This suggests dosage sensitivity in both directions")
    
    # Save results
    gene_effects_df = pd.DataFrame({
        'gene': all_genes,
        'signed_score': all_scores
    })
    gene_effects_df.to_csv(RESULTS_DIR / "revision3b_signed_gene_effects.csv", index=False)
    variant_results.to_csv(RESULTS_DIR / "revision3b_signed_variant_summary.csv", index=False)
    
    print(f"\nResults saved to {RESULTS_DIR}/revision3b_*.csv")
    
    # Create visualization
    create_visualization(all_scores, RESULTS_DIR / "revision3b_signed_distribution.png")
    
    return direction_stats


if __name__ == "__main__":
    main()
