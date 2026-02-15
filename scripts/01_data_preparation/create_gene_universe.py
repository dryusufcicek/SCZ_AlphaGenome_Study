"""
Create Independent Gene Universe

Extracts all protein-coding genes from GENCODE v43 to create an independent
null distribution for empirical Z-score calculation. This ensures we don't
use the complement of the discovery set as background (circular analysis).

Author: SCZ Functional Genomics Pipeline
Date: 2026-02-14
"""

import sys
from pathlib import Path
import pandas as pd
import gzip

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
import config

def parse_gtf_attributes(attr_string):
    """Parse GTF attribute string into dictionary."""
    attrs = {}
    for item in attr_string.strip().split(';'):
        item = item.strip()
        if item:
            key, value = item.split(' ', 1)
            attrs[key] = value.strip('"')
    return attrs

def extract_protein_coding_genes(gtf_file):
    """
    Extract protein-coding genes from GENCODE GTF.

    Parameters
    ----------
    gtf_file : Path
        GENCODE GTF file (gzipped)

    Returns
    -------
    pd.DataFrame
        Gene universe with columns: gene_id, gene_name, chr, start, end, strand
    """
    print("="*80)
    print("EXTRACTING PROTEIN-CODING GENES FROM GENCODE")
    print("="*80)

    print(f"\nParsing GTF file: {gtf_file}")
    genes = []

    with gzip.open(gtf_file, 'rt') as f:
        for line_num, line in enumerate(f, 1):
            # Skip comments
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]

            # Only process gene features
            if feature_type != 'gene':
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = parse_gtf_attributes(fields[8])

            # Only keep protein-coding genes on standard chromosomes
            gene_type = attributes.get('gene_type', attributes.get('gene_biotype', ''))

            if gene_type == 'protein_coding':
                # Only keep autosomal and X chromosome genes
                if chrom.startswith('chr'):
                    chrom = chrom[3:]  # Remove 'chr' prefix

                if chrom in [str(i) for i in range(1, 23)] + ['X']:
                    genes.append({
                        'gene_id': attributes.get('gene_id', ''),
                        'gene_name': attributes.get('gene_name', ''),
                        'chr': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'gene_type': gene_type
                    })

            # Progress
            if line_num % 100000 == 0:
                print(f"  Processed {line_num:,} lines, found {len(genes):,} protein-coding genes...", end='\r')

    print(f"\n  ✓ Parsed {line_num:,} lines total")

    # Create DataFrame
    df = pd.DataFrame(genes)
    print(f"  ✓ Extracted {len(df):,} protein-coding genes")

    # Remove duplicates (keep first occurrence)
    df = df.drop_duplicates(subset=['gene_id'], keep='first')
    print(f"  ✓ {len(df):,} unique genes after deduplication")

    return df

def main():
    """Create gene universe."""
    # Input: GENCODE GTF
    gtf_file = config.DATA_DIR / "external" / "gencode.v43.basic.annotation.gtf.gz"

    if not gtf_file.exists():
        print(f"❌ ERROR: GENCODE GTF file not found: {gtf_file}")
        return 1

    # Extract protein-coding genes
    gene_universe = extract_protein_coding_genes(gtf_file)

    # Summary statistics
    print(f"\nGene Universe Summary:")
    print(f"  Total genes: {len(gene_universe):,}")
    print(f"  Chromosomes: {sorted(gene_universe['chr'].unique())}")
    print(f"\n  Genes per chromosome:")
    chr_counts = gene_universe['chr'].value_counts().sort_index()
    for chr_name, count in chr_counts.items():
        print(f"    chr{chr_name:2s}: {count:5,} genes")

    # Save gene universe
    output_file = config.PROCESSED_DIR / "annotations" / "gene_universe_protein_coding.csv"
    output_file.parent.mkdir(parents=True, exist_ok=True)

    print(f"\nSaving gene universe to: {output_file}")
    gene_universe.to_csv(output_file, index=False)

    print("\n" + "="*80)
    print("✓ GENE UNIVERSE CREATION COMPLETE")
    print("="*80)
    print(f"\nThis gene universe will be used as the independent null")
    print(f"distribution for empirical Z-score calculation.")
    print(f"\nKey properties:")
    print(f"  - Independent of SCZ GWAS discovery")
    print(f"  - All protein-coding genes from GENCODE v43")
    print(f"  - Standard chromosomes only (autosomes + X)")
    print(f"  - {len(gene_universe):,} genes total")

    return 0

if __name__ == "__main__":
    sys.exit(main())
