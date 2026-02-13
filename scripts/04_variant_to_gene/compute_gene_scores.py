#!/usr/bin/env python3
"""Step 4: PP-weighted variant-to-gene aggregation (multi-modal)."""

from __future__ import annotations

import pandas as pd
from pathlib import Path
import sys
from tqdm import tqdm

sys.path.append(str(Path(__file__).parent.parent.parent))
from config import PROCESSED_DIR

MODALITIES = ["DNase", "H3K27ac", "H3K4me1", "H3K4me3", "CAGE", "CTCF", "RNA"]


def resolve_modality_columns(df: pd.DataFrame) -> dict[str, str]:
    """Map canonical modality names to actual score columns in input CSV."""
    lookup = {c.lower(): c for c in df.columns}
    mapping = {}
    for modality in MODALITIES:
        wanted = f"score_{modality}".lower()
        mapping[modality] = lookup.get(wanted)
    return mapping


def aggregate_scores(scores_file: Path) -> pd.DataFrame | None:
    print(f"Loading scores from {scores_file}...")

    try:
        df = pd.read_csv(scores_file)
    except Exception as exc:
        print(f"Failed to read scores file: {exc}")
        return None

    print(f"Loaded {len(df)} variant rows.")

    if "SNP" not in df.columns:
        print("Critical error: SNP column missing in variant score file.")
        return None

    if "PP" not in df.columns:
        print("PP column missing, merging from official credible sets.")
        meta_file = PROCESSED_DIR / "credible_sets" / "scz_credible_sets_official.csv"
        if not meta_file.exists():
            print("Critical error: no PP metadata available.")
            return None
        meta = pd.read_csv(meta_file, usecols=["SNP", "PP", "Locus_ID"])
        df = df.merge(meta, on="SNP", how="left")

    df["PP"] = pd.to_numeric(df["PP"], errors="coerce").fillna(0.0)
    df["target_gene_score"] = pd.to_numeric(df.get("target_gene_score", 0.0), errors="coerce").fillna(0.0)

    mod_cols = resolve_modality_columns(df)
    found_cols = {k: v for k, v in mod_cols.items() if v is not None}
    print(f"Resolved modality columns: {found_cols}")

    gene_data: dict[str, dict[str, float]] = {}

    print("Aggregating PP-weighted gene scores...")
    for _, row in tqdm(df.iterrows(), total=len(df)):
        pp = float(row["PP"])
        if pp <= 0:
            continue

        target_gene = row.get("target_gene")
        if pd.isna(target_gene) or str(target_gene).strip() == "":
            continue
        target_gene = str(target_gene)

        if target_gene not in gene_data:
            gene_data[target_gene] = {
                "weighted_score": 0.0,
                "max_score": 0.0,
                "n_variants": 0,
            }
            for modality in MODALITIES:
                gene_data[target_gene][f"score_{modality}"] = 0.0

        effect_max = float(row.get("target_gene_score", 0.0) or 0.0)
        gene_data[target_gene]["weighted_score"] += effect_max * pp
        gene_data[target_gene]["max_score"] = max(gene_data[target_gene]["max_score"], effect_max)
        gene_data[target_gene]["n_variants"] += 1

        for modality in MODALITIES:
            src_col = mod_cols.get(modality)
            if src_col is None:
                continue
            val = pd.to_numeric(row.get(src_col), errors="coerce")
            if pd.notna(val):
                gene_data[target_gene][f"score_{modality}"] += float(val) * pp

    if not gene_data:
        print("No gene-level entries were generated.")
        return None

    records = []
    for gene, data in gene_data.items():
        rec = {
            "gene": gene,
            "weighted_score": data["weighted_score"],
            "max_score": data["max_score"],
            "n_variants": int(data["n_variants"]),
        }
        for modality in MODALITIES:
            rec[f"score_{modality}"] = data[f"score_{modality}"]
        records.append(rec)

    out_df = pd.DataFrame(records).sort_values("weighted_score", ascending=False)

    print("Top 10 genes by weighted_score:")
    print(out_df[["gene", "weighted_score", "n_variants"]].head(10).to_string(index=False))

    return out_df


def main() -> None:
    print("=" * 70)
    print("STEP 4: PP-WEIGHTED GENE SCORE AGGREGATION")
    print("=" * 70)

    scores_file = PROCESSED_DIR / "alphagenome_variant_scores.csv"
    if not scores_file.exists():
        print(f"Input file missing: {scores_file}")
        return

    res_df = aggregate_scores(scores_file)
    if res_df is None:
        return

    out_file = PROCESSED_DIR / "gene_scores" / "gene_weighted_scores.csv"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    res_df.to_csv(out_file, index=False)
    print(f"Saved {len(res_df)} gene scores to {out_file}")


if __name__ == "__main__":
    main()
