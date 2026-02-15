"""
Posterior Probability Normalization Script

Normalizes posterior probabilities (PP) to sum to 1.0 within each locus.
This fixes the issue where 32.5% of loci have PP sums > 1.0 (up to 5.0),
which violates Bayesian probability axioms and biases gene-level aggregation.

Author: SCZ Functional Genomics Pipeline
Date: 2026-02-14
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
import config

def normalize_posterior_probabilities(input_file, output_file):
    """
    Normalize posterior probabilities within each locus to sum to 1.0.

    Parameters
    ----------
    input_file : Path
        Input credible sets CSV with PP column
    output_file : Path
        Output CSV with PP_normalized column added

    Returns
    -------
    pd.DataFrame
        Normalized credible sets
    """
    print("="*80)
    print("POSTERIOR PROBABILITY NORMALIZATION")
    print("="*80)

    # Load credible sets
    print(f"\nLoading credible sets from: {input_file}")
    df = pd.read_csv(input_file)
    print(f"  Variants: {len(df):,}")
    print(f"  Loci: {df['Locus_ID'].nunique()}")

    # Calculate PP sums per locus BEFORE normalization
    pp_sums_before = df.groupby('Locus_ID')['PP'].sum()
    print(f"\nPP sums BEFORE normalization:")
    print(f"  Min:    {pp_sums_before.min():.6f}")
    print(f"  Median: {pp_sums_before.median():.6f}")
    print(f"  Mean:   {pp_sums_before.mean():.6f}")
    print(f"  Max:    {pp_sums_before.max():.6f}")

    problem_loci = (pp_sums_before > 1.0).sum()
    print(f"  Loci with PP > 1.0: {problem_loci} / {len(pp_sums_before)} ({100*problem_loci/len(pp_sums_before):.1f}%)")

    # Identify worst offenders
    worst_loci = pp_sums_before.nlargest(5)
    print(f"\n  Top 5 problematic loci:")
    for locus_id, pp_sum in worst_loci.items():
        n_variants = (df['Locus_ID'] == locus_id).sum()
        locus_name = df.loc[df['Locus_ID'] == locus_id, 'Locus_Name'].iloc[0]
        print(f"    Locus {locus_id:3d} ({locus_name:15s}): PP sum = {pp_sum:.4f} ({n_variants} variants)")

    # Normalize PP within each locus
    print(f"\nNormalizing PP within each locus...")
    df['PP_normalized'] = df.groupby('Locus_ID')['PP'].transform(lambda x: x / x.sum())

    # Verify normalization
    pp_sums_after = df.groupby('Locus_ID')['PP_normalized'].sum()
    print(f"\nPP sums AFTER normalization:")
    print(f"  Min:    {pp_sums_after.min():.6f}")
    print(f"  Median: {pp_sums_after.median():.6f}")
    print(f"  Mean:   {pp_sums_after.mean():.6f}")
    print(f"  Max:    {pp_sums_after.max():.6f}")

    # Check for any deviations from 1.0 (should be none, or only floating point errors)
    deviations = np.abs(pp_sums_after - 1.0)
    max_deviation = deviations.max()
    print(f"  Max deviation from 1.0: {max_deviation:.2e}")

    if max_deviation > 1e-10:
        print(f"  ⚠️  WARNING: Some loci have PP sum != 1.0 (max error: {max_deviation:.2e})")
    else:
        print(f"  ✓ All loci normalized correctly (max floating-point error: {max_deviation:.2e})")

    # Keep original PP for reference
    df = df.rename(columns={'PP': 'PP_original'})
    df = df.rename(columns={'PP_normalized': 'PP'})

    # Save normalized credible sets
    print(f"\nSaving normalized credible sets to: {output_file}")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, index=False)
    print(f"  ✓ Saved {len(df):,} variants")

    # Generate comparison statistics
    print(f"\nNormalization Impact:")
    # For loci that had PP > 1.0, how much were individual variant PPs reduced?
    problem_locus_ids = pp_sums_before[pp_sums_before > 1.0].index
    problem_variants = df[df['Locus_ID'].isin(problem_locus_ids)]

    if len(problem_variants) > 0:
        reduction_factors = problem_variants['PP_original'] / problem_variants['PP']
        print(f"  Affected variants: {len(problem_variants):,} ({100*len(problem_variants)/len(df):.1f}%)")
        print(f"  PP reduction factors:")
        print(f"    Min (least reduced):  {reduction_factors.min():.4f}x")
        print(f"    Median:               {reduction_factors.median():.4f}x")
        print(f"    Max (most reduced):   {reduction_factors.max():.4f}x")

    return df

def main():
    """Run PP normalization."""
    # Input: original credible sets with unnormalized PP
    input_file = config.CS_DIR / "scz_credsets_ST11d_grch38_pip.csv"

    # Output: normalized credible sets (PP now sums to 1.0 per locus)
    output_file = config.CS_DIR / "scz_credsets_normalized.csv"

    # Normalize
    df = normalize_posterior_probabilities(input_file, output_file)

    print("\n" + "="*80)
    print("✓ POSTERIOR PROBABILITY NORMALIZATION COMPLETE")
    print("="*80)
    print(f"\nNormalized credible sets saved to:")
    print(f"  {output_file}")
    print(f"\nUse this file for all downstream analyses to ensure")
    print(f"equal weighting across loci in gene aggregation.")

    return 0

if __name__ == "__main__":
    sys.exit(main())
