"""
Extract Brain eGenes from GTEx v10

Extracts genes with significant brain eQTLs from GTEx v10 eGenes files.
Creates a comprehensive list of genes with brain eQTL evidence for prioritization.

Author: SCZ Functional Genomics Pipeline
Date: 2026-02-15
"""

import sys
from pathlib import Path
import pandas as pd
import gzip

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
import config

# GTEx directory
GTEX_DIR = config.PROCESSED_DIR / "eqtl_data" / "GTEx_Analysis_v10_eQTL_updated"

# 13 brain tissues in GTEx v10
BRAIN_TISSUES = [
    'Brain_Amygdala',
    'Brain_Anterior_cingulate_cortex_BA24',
    'Brain_Caudate_basal_ganglia',
    'Brain_Cerebellar_Hemisphere',
    'Brain_Cerebellum',
    'Brain_Cortex',
    'Brain_Frontal_Cortex_BA9',
    'Brain_Hippocampus',
    'Brain_Hypothalamus',
    'Brain_Nucleus_accumbens_basal_ganglia',
    'Brain_Putamen_basal_ganglia',
    'Brain_Spinal_cord_cervical_c-1',
    'Brain_Substantia_nigra',
]

def load_egenes_from_tissue(tissue_name):
    """
    Load eGenes from a single brain tissue.

    Parameters
    ----------
    tissue_name : str
        GTEx tissue name (e.g., 'Brain_Cortex')

    Returns
    -------
    pd.DataFrame
        eGenes with columns: gene_id, gene_name, tissue, qval
    """
    egenes_file = GTEX_DIR / f"{tissue_name}.v10.eGenes.txt.gz"

    if not egenes_file.exists():
        print(f"⚠️  {tissue_name} eGenes file not found")
        return None

    # Read gzipped file
    df = pd.read_csv(egenes_file, sep='\t', compression='gzip')

    # Filter for significant eGenes (qval < 0.05)
    egenes = df[df['qval'] < 0.05].copy()

    # Keep relevant columns
    egenes = egenes[['gene_id', 'gene_name', 'qval', 'pval_nominal']].copy()
    egenes['tissue'] = tissue_name

    return egenes

def extract_all_brain_egenes():
    """
    Extract eGenes from all 13 brain tissues.

    Returns
    -------
    pd.DataFrame
        Combined eGenes from all brain tissues
    """
    print("="*80)
    print("EXTRACTING BRAIN eGENES FROM GTEx v10")
    print("="*80)

    print(f"\nGTEx directory: {GTEX_DIR}")
    print(f"Brain tissues: {len(BRAIN_TISSUES)}")

    all_egenes = []

    for tissue in BRAIN_TISSUES:
        print(f"\nProcessing {tissue}...")

        egenes = load_egenes_from_tissue(tissue)

        if egenes is None:
            continue

        n_egenes = len(egenes)
        print(f"  ✓ {n_egenes:,} eGenes (qval < 0.05)")

        all_egenes.append(egenes)

    # Combine all tissues
    combined_egenes = pd.concat(all_egenes, ignore_index=True)

    print(f"\n{'='*80}")
    print(f"COMBINED eGENES ACROSS ALL BRAIN TISSUES")
    print(f"{'='*80}")

    print(f"\nTotal eGene records: {len(combined_egenes):,}")
    print(f"Unique genes: {combined_egenes['gene_name'].nunique():,}")

    # Count tissues per gene
    tissues_per_gene = combined_egenes.groupby('gene_name')['tissue'].nunique()

    print(f"\nTissue coverage per gene:")
    print(f"  Min: {tissues_per_gene.min()} tissues")
    print(f"  Max: {tissues_per_gene.max()} tissues")
    print(f"  Mean: {tissues_per_gene.mean():.1f} tissues")
    print(f"  Median: {tissues_per_gene.median():.1f} tissues")

    # Genes with eQTLs in many tissues
    highly_supported = tissues_per_gene[tissues_per_gene >= 5]
    print(f"\nGenes with eQTLs in ≥5 brain tissues: {len(highly_supported):,}")

    return combined_egenes

def create_gene_eqtl_lookup(combined_egenes):
    """
    Create gene-level eQTL lookup table.

    For each gene, summarize:
    - Whether it has brain eQTL evidence
    - Number of brain tissues with eQTL
    - Minimum qval across tissues

    Parameters
    ----------
    combined_egenes : pd.DataFrame
        Combined eGenes from all tissues

    Returns
    -------
    pd.DataFrame
        Gene-level eQTL summary
    """
    print("\n" + "="*80)
    print("CREATING GENE-LEVEL eQTL LOOKUP")
    print("="*80)

    # Aggregate by gene
    gene_summary = combined_egenes.groupby('gene_name').agg({
        'tissue': 'count',  # Number of tissues with eQTL
        'qval': 'min',      # Best (minimum) q-value
    }).reset_index()

    gene_summary.columns = ['gene_symbol', 'n_tissues_with_eqtl', 'min_qval']

    # Add flag
    gene_summary['has_brain_eqtl'] = True

    print(f"\nGene-level summary:")
    print(f"  Total genes with brain eQTL: {len(gene_summary):,}")

    # Distribution
    print(f"\nDistribution of tissue coverage:")
    tissue_dist = gene_summary['n_tissues_with_eqtl'].value_counts().sort_index()
    for n_tissues, count in tissue_dist.items():
        print(f"  {n_tissues:2d} tissues: {count:5,} genes")

    # Save
    output_file = config.PROCESSED_DIR / "eqtl_data" / "gtex_brain_egenes_summary.csv"
    gene_summary.to_csv(output_file, index=False)

    print(f"\n✓ Saved gene-level eQTL lookup: {output_file.name}")

    # Also save full detailed table
    detail_file = config.PROCESSED_DIR / "eqtl_data" / "gtex_brain_egenes_detailed.csv"
    combined_egenes.to_csv(detail_file, index=False)

    print(f"✓ Saved detailed eGenes table: {detail_file.name}")

    return gene_summary

def check_known_scz_genes(gene_summary):
    """Check eQTL coverage for known SCZ genes."""
    print("\n" + "="*80)
    print("eQTL COVERAGE FOR KNOWN SCZ GENES")
    print("="*80)

    known_scz_genes = [
        'CACNA1C', 'TCF4', 'ZNF804A', 'NRGN', 'GRIN2A', 'GRM3',
        'DRD2', 'SRR', 'CLCN3', 'SNAP91', 'FURIN', 'TSNARE1',
        'AS3MT', 'CNNM2', 'NT5C2', 'ITIH3', 'CNTN4', 'PCCB'
    ]

    print(f"\nChecking {len(known_scz_genes)} known SCZ genes:")
    print(f"{'Gene':12s} {'Brain eQTL?':>12s} {'Tissues':>8s} {'Min qval':>10s}")
    print("-" * 50)

    for gene in known_scz_genes:
        gene_data = gene_summary[gene_summary['gene_symbol'] == gene]

        if len(gene_data) > 0:
            n_tissues = gene_data['n_tissues_with_eqtl'].values[0]
            min_qval = gene_data['min_qval'].values[0]
            print(f"{gene:12s} {'✓':>12s} {n_tissues:8d} {min_qval:10.2e}")
        else:
            print(f"{gene:12s} {'✗':>12s} {'N/A':>8s} {'N/A':>10s}")

    # Summary
    n_with_eqtl = sum(1 for gene in known_scz_genes if gene in gene_summary['gene_symbol'].values)
    print(f"\n{n_with_eqtl}/{len(known_scz_genes)} known SCZ genes have brain eQTL evidence ({100*n_with_eqtl/len(known_scz_genes):.1f}%)")

def main():
    """Main execution."""
    print("╔" + "="*78 + "╗")
    print("║" + " "*20 + "GTEx v10 BRAIN eGENES EXTRACTION" + " "*26 + "║")
    print("╚" + "="*78 + "╝")

    # Step 1: Extract eGenes from all brain tissues
    combined_egenes = extract_all_brain_egenes()

    # Step 2: Create gene-level lookup
    gene_summary = create_gene_eqtl_lookup(combined_egenes)

    # Step 3: Check known SCZ genes
    check_known_scz_genes(gene_summary)

    print("\n" + "="*80)
    print("✓ GTEx BRAIN eGENES EXTRACTION COMPLETE")
    print("="*80)

    print(f"""
SUMMARY:
- Extracted eGenes from 13 brain tissues
- Total unique genes with brain eQTL: {len(gene_summary):,}
- Known SCZ genes with brain eQTL: see table above

FILES CREATED:
- {config.PROCESSED_DIR / "eqtl_data" / "gtex_brain_egenes_summary.csv"}
  (gene-level: gene_symbol, n_tissues, min_qval, has_brain_eqtl)

- {config.PROCESSED_DIR / "eqtl_data" / "gtex_brain_egenes_detailed.csv"}
  (tissue-level: all eGene-tissue pairs)

NEXT STEP:
Use gtex_brain_egenes_summary.csv to prioritize genes in aggregation.
    """)

    return 0

if __name__ == "__main__":
    sys.exit(main())
