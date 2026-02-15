# Main Figure Legends

## Figure 1. Multi-Layer Regulatory Architecture of Schizophrenia Risk Genes

**A. Study design overview.** Integration pipeline combining fine-mapped PGC3 schizophrenia GWAS variants (20,760 variants, 255 loci) with AlphaGenome deep learning regulatory predictions (six chromatin modalities) and brain-specific 3D chromatin architecture (PsychENCODE fetal Hi-C). Gene-level aggregation incorporates both linear genomic windows (±256 kb) and long-range chromatin loops (up to 1 Mb). Statistical enrichment testing against independent gene universe (19,966 protein-coding genes) identifies 1,617 genes with significant regulatory impact (FDR<0.10). Independent validation via GTEx brain eQTL (13 tissues).

**B. Regulatory gene discovery by mapping source.** Venn diagram showing genes identified through linear window mapping (3,725 genes), Hi-C chromatin loops (502 genes), and their overlap (283 genes with both sources). 219 genes (5.6% of total) are Hi-C-exclusive, accessible only through distal chromatin contacts. Bar chart (right) shows eQTL validation rates: Hi-C genes show 84% brain eQTL concordance vs 44.6% for linear-only genes (1.88-fold enrichment, Fisher's exact p<0.0001).

**C. Brain-region specificity of regulatory enrichment.** Heatmap showing eQTL enrichment across 13 GTEx brain tissues (rows) for ranked gene tiers (columns: top 100, 101–500, 501–1,000, 1,001–1,617). Color intensity: -log₁₀(p-value) from Mann-Whitney U test. Strongest enrichment in frontal cortex (BA9, p=0.000366), cerebellum (p=0.000465), and nucleus accumbens (p=0.000701). 12 of 13 tissues show significant enrichment (FDR<0.05).

**D. Modality-specific regulatory signatures.** Radar plots showing normalized z-scores across six chromatin modalities for representative gene archetypes: **Chromatin-dominant** (CCNT2: high DNase/ATAC, low RNA), **Transcription initiation-dominant** (TBC1D17: high CAGE/PROCAP), **RNA-dominant** (NAPRT: high RNA, pan-brain eQTL), **Splicing-specific** (KLF16: high splice junction). eQTL concordance varies by modality: RNA-dominant (75%), chromatin-dominant (60%), splicing-specific (45%).

---

## Figure 2. 3D Chromatin Architecture Reveals Distal Regulatory Connections

**A. Hi-C loop integration overview.** Schematic showing variant-to-gene mapping via chromatin loops. Top: Linear genomic view showing schizophrenia variant (red star) within ±256 kb window captures three proximal genes (blue). Bottom: 3D chromatin architecture view showing same variant in loop anchor contacts distal gene promoter 600 kb away (green), accessible only through spatial genome folding. PsychENCODE fetal brain Hi-C: 149,097 loops tested, 3,269 variants in loop anchors.

**B. Distance distribution for Hi-C connections.** Histogram showing genomic distances for 5,941 Hi-C variant-gene pairs. Median: 556 kb, range 10–1,019 kb. Shaded region (±256 kb): only 15% of Hi-C connections fall within linear window range. 85% of Hi-C connections are >256 kb, validating recovery of long-range regulatory relationships.

**C. Top Hi-C-exclusive genes.** Bar chart showing composite scores for top 20 Hi-C-exclusive genes (accessible only through chromatin loops, no linear window connection). Top genes: ZNF554 (rank 1, score 5.21, distance >700 kb), KLF16 (rank 3, score 2.24), S100A8 (rank 4, score 1.63, immune response), TBC1D17 (rank 9, score 0.97, vesicular trafficking), CCNT2 (rank 11, score 0.92, transcriptional elongation). Color intensity: number of Hi-C variants connecting to gene.

**D. Multi-scale regulatory architecture.** Scatter plot showing genes with both linear and Hi-C evidence (n=283 genes). X-axis: number of linear variants. Y-axis: number of Hi-C variants. Point size: composite score. Representative genes labeled: SEZ6L2 (36 linear + 24 Hi-C, synaptic adhesion), KCTD13 (36 linear + 9 Hi-C, 16p11.2 deletion syndrome), MVP (36 linear + 18 Hi-C, nuclear transport). Genes with convergent evidence from multiple spatial scales may represent critical regulatory nodes.

---

## Figure 3. Brain-Region Specificity and Regulatory Validation

**A. Tissue-specific eQTL enrichment.** Forest plot showing fold enrichment for brain eQTL in top 1,617 genes versus background (all 19,966 genes, 48.8% eQTL rate). Tissues ordered by effect size. Frontal cortex (BA9): 1.82-fold enrichment (95% CI 1.45–2.28, p=0.000366). Cerebellum: 1.76-fold (1.41–2.19, p=0.000465). Nucleus accumbens: 1.71-fold (1.38–2.12, p=0.000701). 12 of 13 tissues FDR-significant (dashed line: p=0.05).

**B. Pan-brain vs region-restricted regulatory programs.** Left: Histogram showing distribution of tissue breadth (number of tissues with eQTL per gene). 367 genes (23%) show pan-brain regulation (≥10 tissues). Representative pan-brain genes: CCNT2 (13/13 tissues), NAPRT (13/13), ATF5 (12/13). Right: Region-restricted genes (1–2 tissues, n=20 in top 100): ACMSD (spinal cord only, tryptophan metabolism), NR1H2 (substantia nigra only, lipid metabolism), NAPSA (anterior cingulate only).

**C. eQTL validation by gene ranking.** Line plot showing percentage of genes with brain eQTL (Y-axis) across ranking tiers (X-axis): top 100 (69%, 2.33-fold enrichment vs background, p<0.0001), 101–500 (61%, 1.47-fold, p=0.0003), 501–1,000 (54%, 1.25-fold, p=0.021), 1,001–1,617 (48%, 0.98-fold, p=0.52). Shaded region: 95% CI. Dashed line: genome-wide background rate (48.8%). Progressive enrichment validates that composite scores capture real regulatory effects.

**D. Known schizophrenia genes: validation and coverage.** Bar chart showing ranks for established schizophrenia genes. Successfully prioritized: CACNA1C (rank 1,455, FDR=0.081, marginally significant, 82 linear variants), CACNB2 (rank 1,345, improved by Hi-C: 54 linear + 2 Hi-C variants), GRIN2A (rank 2,031, weak signal). Null results (fetal Hi-C limitation): DRD2 (rank 2,350, 0 Hi-C variants), GRM3 (rank 2,796, 0 Hi-C variants), GAD1/GAD2/SYN1/SYN2 (not in dataset). Bar color indicates mapping source (blue: linear only, green: linear+Hi-C, red: null result).

---

## Figure 4. Pathway-Level Organization and Novel Gene Discovery

**A. Gene Set Enrichment Analysis results.** Dot plot showing top 20 enriched pathways (nominal p<0.05). X-axis: mean rank of pathway genes (lower = higher scores). Y-axis: pathway name. Point size: number of genes. Color: -log₁₀(p-value). Top synaptic pathways: regulation of calcium ion-dependent exocytosis (p=0.0007, 2 genes), regulation of short-term synaptic plasticity (p=0.0043, 1 gene), regulation of synapse organization (p=0.0224, 2 genes). Top non-synaptic: defense response to bacterium (p=0.000048, 5 genes, immune/inflammation), amino acid catabolism (p=0.00024, 4 genes). Dashed line: FDR q=0.10 threshold. **None FDR-significant.**

**B. Missing synaptic genes: fetal Hi-C coverage.** Schematic showing canonical schizophrenia synaptic genes absent from Hi-C dataset. Top: DRD2 locus showing 69 linear variants (blue dots) within ±256 kb but no Hi-C loop anchors (gray box: no loops in this region). Bottom: PsychENCODE fetal Hi-C loop coverage across genome. Red regions: high loop density. White gaps: no loops. Canonical synaptic gene loci (DRD2, GRM3, GAD1, GAD2) fall in white gaps, explaining null results. Inset: temporal specificity hypothesis—fetal brain Hi-C captures early neurodevelopmental programs; adult synaptic genes regulated through later-established architectures.

**C. Top novel schizophrenia risk gene candidates.** Heatmap showing top 20 novel genes (high rank, no prior strong schizophrenia association) × six regulatory modalities. Rows: genes ordered by rank. Columns: modality z-scores (color: blue=low, red=high). Annotations (right): brain eQTL tissues, mapping source (linear/Hi-C), gene function. Top candidates: ZNF554 (rank 1, Hi-C-exclusive, transcriptional regulator), KLF16 (rank 3, splicing-dominant, transcription factor), S100A8 (rank 4, Hi-C-exclusive, immune/calcium signaling), TBC1D17 (rank 9, transcription initiation-dominant, vesicular trafficking), CCNT2 (rank 11, chromatin-dominant, pan-brain eQTL, transcriptional elongation).

**D. Multi-layer regulatory framework summary.** Conceptual diagram integrating key findings. Center: schizophrenia risk variant (red). Layer 1: Local chromatin architecture (±256 kb, 3,725 genes, 44.6% eQTL validation). Layer 2: Distal chromatin loops (up to 1 Mb, 502 genes including 219 Hi-C-exclusive, 84% eQTL validation). Layer 3: Modality-specific regulation (chromatin accessibility, transcription initiation, RNA abundance, splicing). Layer 4: Brain region specificity (strongest: frontal cortex, cerebellum, nucleus accumbens; moderate: hippocampus, amygdala; pan-brain: 367 genes). Outer ring: therapeutic implications (chromatin modifiers, transcription factors, RNA stability regulators, splicing modulators).

---

## Extended Data Figure Legends

### Extended Data Figure 1. Fine-Mapping Quality Control

**A-D.** Posterior probability distributions and normalization validation. **A:** Per-locus PP sum distribution (n=255 loci). **B:** Before-after normalization histograms. **C:** Variant-level PP changes (scatter plot, n=20,760 variants). **D:** Gene rank impact of PP normalization (top 1,000 genes, color-coded by rank change magnitude).

---

### Extended Data Figure 2. AlphaGenome Prediction Characteristics

**A-D.** Modality-specific score properties. **A:** Raw score distributions (violin plots, log scale). **B:** Z-scored distributions (violin plots). **C:** Correlation heatmap (6×6 modalities). **D:** Example gene profiles showing modality-specific signatures (n=12 representative genes).

---

### Extended Data Figure 3. Hi-C Integration Technical Validation

**A-D.** Hi-C mapping validation and characteristics. **A:** Loop anchor size distribution. **B:** Variant-anchor overlap statistics. **C:** H-MAGMA concordance analysis (93.2% agreement). **D:** Hi-C gene enrichment in top-ranked genes (Fisher's exact tests for top 100, 500, 1,000).

---

### Extended Data Figure 4. Brain eQTL Validation Details

**A-D.** Comprehensive eQTL validation. **A:** eQTL coverage per tissue (sample sizes and eGene counts). **B:** Tissue-specific enrichment heatmaps (13 tissues × rank tiers). **C:** eQTL effect size distributions (discovery genes vs background). **D:** Genes lacking eQTL: characteristics and potential explanations (cell type-specific, developmental, post-transcriptional).

---

### Extended Data Figure 5. Pathway Enrichment Comprehensive Results

**A-D.** Full pathway analysis. **A:** Volcano plot (all 3,623 pathways tested). **B:** Top 50 pathways by nominal p-value. **C:** Custom neuroscience pathways (all p>0.05). **D:** Gene overlap between top enriched pathways.

---

### Extended Data Figure 6. Known Schizophrenia Genes Detailed Analysis

**A-D.** Established genes: successes and limitations. **A:** Ranks and scores for 25 known schizophrenia genes. **B:** Variant counts (linear vs Hi-C) for known genes. **C:** Modality profiles for successfully prioritized genes (CACNA1C, CACNB2, GRIN2A). **D:** Hi-C coverage analysis explaining null results (DRD2, GRM3, GAD1/GAD2).

---

### Extended Data Figure 7. Novel Gene Functional Characterization

**A-D.** Top novel candidates: evidence summary. **A:** Protein-protein interaction networks for top 20 novel genes. **B:** Gene Ontology enrichment for novel gene set. **C:** Expression patterns across brain development (BrainSpan data). **D:** Single-cell expression profiles in human cortex (excitatory neurons, inhibitory neurons, glia).

---

## Figure File Specifications

**Main Figures 1-4:**
- Format: PDF (vector) with 300 DPI embedded rasters for heatmaps
- Dimensions: Figure 1 (full page, 183×247 mm), Figures 2-4 (half page, 183×120 mm)
- Fonts: Arial or Helvetica, minimum 6pt, labels 8-10pt
- Color scheme: Colorblind-friendly (viridis, RdBu, Set2)
- Line weights: Minimum 0.5pt for visibility at journal print scale

**Extended Data Figures 1-7:**
- Format: PDF (vector) with 300 DPI embedded rasters
- Dimensions: Standard A4 or Letter
- Same font and color specifications as main figures

**Supplementary Figure Files:**
- Format: High-resolution PNG (600 DPI) or PDF
- Deposited with supplementary information

---

## Data Visualization Guidelines

All figures follow Nature Neuroscience data visualization standards:
- **Error bars:** 95% confidence intervals (Fisher's exact tests) or standard error of mean (continuous data)
- **Statistical annotations:** Exact p-values when p<0.001; otherwise p-value ranges (*p<0.05, **p<0.01, ***p<0.001)
- **Multiple testing:** All p-values FDR-corrected where applicable; correction method specified in legends
- **Color accessibility:** All heatmaps and scatter plots use colorblind-safe palettes
- **Sample sizes:** Displayed in figure legends or on plots
- **Reproducibility:** All figures generated via documented scripts (Code Availability section)

---

## Figure Source Data

**Figure 1:**
- Source Data Fig. 1B: Gene counts by mapping source, eQTL validation rates
- Source Data Fig. 1C: eQTL enrichment p-values (13 tissues × 4 rank tiers)
- Source Data Fig. 1D: Modality z-scores for representative genes (n=4)

**Figure 2:**
- Source Data Fig. 2B: Hi-C distance distribution (n=5,941 variant-gene pairs)
- Source Data Fig. 2C: Top 20 Hi-C-exclusive genes with scores and distances
- Source Data Fig. 2D: Multi-scale genes (linear + Hi-C variant counts, n=283)

**Figure 3:**
- Source Data Fig. 3A: Tissue-specific fold enrichment and 95% CI (n=13)
- Source Data Fig. 3B: Pan-brain gene examples, region-restricted gene examples
- Source Data Fig. 3C: eQTL rates by rank tier (4 tiers)
- Source Data Fig. 3D: Known gene ranks and variant counts (n=25 genes)

**Figure 4:**
- Source Data Fig. 4A: Pathway enrichment statistics (top 20 pathways)
- Source Data Fig. 4B: Hi-C loop coverage at canonical synaptic gene loci
- Source Data Fig. 4C: Top 20 novel genes × 6 modalities (z-score matrix)
- Source Data Fig. 4D: Multi-layer framework summary statistics

**All source data deposited in Supplementary Tables S1-S7.**

