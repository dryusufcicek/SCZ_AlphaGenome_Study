# Supplementary Information

**Title:** Deep Learning Reveals the Regulatory Architecture of Schizophrenia Risk
**Authors:** Yusuf Cicek, et al.

---

## Supplementary Note 1: GWAS Data and Locus Definition
**Source:** We utilized the Psychiatric Genomics Consortium Wave 3 (PGC3) schizophrenia summary statistics (Trubetskoy et al., 2022), comprising 76,755 cases and 243,649 controls of predominantly European ancestry.
**Locus Definition:** Genomic loci were defined by iteratively merging genome-wide significant variants ($P < 5 \times 10^{-8}$) located within 500kb of each other. This procedure yielded **237 independent genomic loci**.
**Filtering:** Variants in the Major Histocompatibility Complex (MHC) region (chr6:25-34Mb) were excluded from primary scoring due to complex LD structure, though they were retained for specific sensitivity analyses.
**Reference Panel:** Linkage Disequilibrium (LD) calculations for credible set construction utilized the 1000 Genomes Project Phase 3 reference panel (European super-population).

## Supplementary Note 2: Credible Set Construction
**Methodology:** To prioritize causal variants within each locus, we calculated Approximate Bayes Factors (ABF) assuming a single causal variant model per block as a simplifying heuristic for large-scale scoring.
**Posterior Probabilities (PP):** For each variant $i$ in locus $L$, the posterior probability was calculated as:
$$PP_i = \frac{ABF_i}{\sum_{j \in L} ABF_j}$$
**Credible Sets:** We constructed 95% Credible Sets by ranking variants by $PP$ and retaining the cumulative sum up to 0.95.
**Summary Statistics:**
*   Total prioritized variants: **10,832**
*   Mean variants per credible set: **45.7**
*   Median variants per credible set: **22**
This approach ensures that our analysis captures the full regulatory potential of the locus, rather than relying on a single, potentially non-causal lead SNP.

## Supplementary Note 3: AlphaGenome Model Details
**Model Architecture:** AlphaGenome is a sequence-based deep neural network trained to predict epigenetic profiles from DNA sequence. We utilized the standard pre-trained checkpoint available via the DeepMind API.
**Input/Output:**
*   **Input context:** 131,072 bp (128kb) centered on the variant.
*   **Output modalities:** DNAse-seq, ATAC-seq, H3K27ac, H3K4me3, H3K4me1.
*   **Biosamples:** Predictions were generated for **53** relevant biosamples, including adult brain regions (Cortex, Hippocampus, Substantia Nigra), fetal brain, and non-neural controls (Liver, Immune, Muscle).
**Disclaimer:** AlphaGenome scores were used strictly as variant-level annotations to inform statistical fine-mapping, not as independent assertions of causality.

## Supplementary Note 4: Variant-to-Gene Scoring Formula
**Aggregation Logic:** To calculate the regulatory burden for a specific gene $g$, we summed the contributions of all variants $v$ in its regulatory vicinity (promoter + enhancers), weighted by their causal probability ($PP_v$).
**Formula:**
$$GeneScore_g = \sum_{v \in CredibleSets} (PP_v \times |LFC_{v,g}|)$$
Where:
*   $PP_v$ is the Posterior Probability of variant $v$ being causal.
*   $LFC_{v,g}$ is the AlphaGenome-predicted Log-Fold Change in chromatin accessibility/activity at gene $g$ caused by variant $v$.
*   We utilized independent summation to capture additive regulatory effects across multiple LD blocks affecting the same gene.

## Supplementary Note 5: Empirical Null & Z-Score Calibration
**Rationale:** Raw deep learning scores are non-Gaussian and difficult to interpret (e.g., is a score of 5.0 high?). We calibrated these against a genome-wide background.
**Background Universe:** We defined a "Null Universe" of **47,808 genes**, representing the comprehensive set of protein-coding and lncRNA genes in the genome.
**Z-Score Formulation:**
$$Z_g = \frac{Score_g - \mu_{null}}{\sigma_{null}}$$
Where $\mu_{null}$ and $\sigma_{null}$ are the mean and standard deviation of scores across the entire 47k gene universe.
**Result:** A Z-score of 4.0 implies the gene carries a regulatory burden 4 standard deviations above the genome-wide expectation. This empirically calibrated metric allows for fair comparison between genes.

## Supplementary Note 6: Sensitivity Analyses
**Robustness Check:** To ensure results were not driven by extreme outliers (Table 1), we performed a "Leave-One-Percent-Out" sensitivity analysis.
**Procedure:**
1.  Identified the top 1% of genes by Z-score ($N \approx 470$).
2.  Removed these genes from the ranked list.
3.  Re-ran GSEA on the remaining 99%.
**Results:** As detailed in the main text, while the hyper-significant P-values for "Transcription" attenuated, the core modules of **Synapse Organization** and **Heterochromatin Organization** remained statistically significant ($P < 0.05$), confirming the polygenic nature of the signal.

## Supplementary Note 7: Unbiased GSEA Methodology
**Pathway Database:** Gene Ontology (GO) Biological Process 2023 (~5,000 terms).
**Statistical Test:** Mann-Whitney U test (one-sided, testing for higher rank).
**Ranking Metric:** Empirical Z-Score. Non-overlapping genes were assigned the minimum rank (0).
**Multiple Testing:** Benjamini-Hochberg False Discovery Rate (FDR) correction was applied across all 5,000 tested pathways.
**Transparency:** The full ranked list of all tested pathways is provided in **Supplementary Table S3**.

## Supplementary Note 8: Cell-Type Enrichment & Footprint Normalization
**Data Source:** Cell-type specific ATAC-seq peaks from the human brain (Corces et al., 2020 / Trevino et al., 2021).
**The "Volume" Problem:** iPSCs and progenitors have significantly larger Total Accessible Chromatin volumes than mature neurons, leading to identifying more variant overlaps by chance.
**Exact Normalization:** We calculated the binomial probability of overlap $P(k; n, p)$, where $p$ was adjusted for total accessible base pairs:
$$p_{celltype} = \frac{\text{Total ATAC bp in Celltype}}{\text{Total Genome Size}}$$
**Comparison:**
*   **Raw Overlap:** iPSC enrichment $P < 10^{-10}$ (Spurious).
*   **Normalized:** iPSC enrichment $P = 0.07$ (Not Significant).
*   **Normalized:** Excitatory Neuron enrichment $P < 10^{-40}$ (Robust).

## Supplementary Note 9: Directionality
**Observation:** Our analysis quantified the *magnitude* of regulatory disruption (|LFC|).
**Global Trend:** Preliminary analysis of signed LFCs suggests a mixed landscape, with no single global skew towards upregulation or downregulation across all risk loci. This is consistent with complex traits where risk can be conferred by both gain-of-function and loss-of-function mechanisms, or by stabilizing/destabilizing chromatin state in context-dependent ways. Thus, we focus on "Dysregulation" magnitude rather than simple directionality.

## Supplementary Note 10: Benchmarking
**Comparison:** We compared AlphaGenome-prioritized genes ($Z > 2$) with:
1.  **Nearest Gene:** Simple proximity mapping.
2.  **PGC3 Prioritized:** Official study list (Trubetskoy et al.).
**Overlap:** We observed a significant overlap (OR = 5.03, $P < 10^{-10}$) with the PGC3 prioritized list, confirming that our model recovers known biology while expanding the candidate set to include mechanistically linked partners (e.g., *H2AC20*) that were missed by statistical fine-mapping alone.

## Supplementary Note 11: Limitations & Model Assumptions
**Adult Tissue Bias:** The AlphaGenome model was primarily trained on adult tissue data (ENCODE/GTEx). While we include fetal brain biosamples, the model may under-detect strictly developmental regulatory events.
**Probabilistic Nature:** All predictions are computational inferences. While statistically calibrated, they have not been validated by direct CRISPR perturbation in this specific study.
**Missing Heritability:** Our credible sets explain the *common variant* signal; rare variants provided by sequencing studies (SCHEMA) offer complementary but distinct insights.

---

## Supplementary Data Tables
**Table S1:** **All Scored Variants.** CSV containing RSID, Chromosome, Position, PGC3 P-value, Proxy Posterior Probability, and Raw AlphaGenome Scores.
**Table S2:** **Gene-Level Results.** CSV containing Gene Symbol, Aggregated Weighted Score, Empirical Z-Score, and Rank.
**Table S3:** **Full GSEA Results.** CSV containing 5,000 GO Biological Processes, N_Overlap, P-value, FDR, and Effect Size.
**Table S4:** **Cell-Type Specificity.** CSV containing raw counts and normalized binomial statistics for all tested cell types.
**Table S5:** **Sensitivity Analysis.** Comparison of top pathways before and after removing top 1% outlier genes.
