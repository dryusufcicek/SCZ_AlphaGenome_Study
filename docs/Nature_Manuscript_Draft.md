# Deep Learning-Based Functional Annotation of Regulatory Risk in Schizophrenia

**Authors:** Yusuf Cicek, et al.
**Affiliation:** Istanbul University - Cerrahpasa

## Abstract
Schizophrenia (SCZ) polygenic risk remains poorly understood at the mechanistic level. Here, we applied **AlphaGenome**, a sequence-based deep learning framework, to decode the regulatory impact of 10,832 candidate causal variants prioritized from the PGC3 GWAS. By implementing a rigorous **Posterior-Weighted Credible Set** analysis and **Unbiased Whole-Genome Enrichment** strategy, we demonstrate that SCZ risk converges fundamentally on the disruption of **neuronal transcriptional regulation** ($P < 10^{-30}$), which is accompanied by specific downstream deficits in **synaptic vesicle transport** ($P < 10^{-24}$) and **intracellular calcium homeostasis** ($P < 10^{-5}$). Crucially, we show via chromatin-footprint normalization that these effects are **predominantly enriched in mature excitatory neurons with additional contributions from glial regulatory programs**, arguing against progenitor-dominant models of risk variance. These findings support a coherent regulatory framework linking non-coding risk to transcriptional and synaptic pathways.

---

## 1. Introduction
Schizophrenia (SCZ) is highly heritable, yet the functional logic of its non-coding risk variants remains obscure. Current methods relying on proximity or broad chromatin interactions often yield non-specific results. Unlike proximity-based or correlation-driven methods, sequence-based models like AlphaGenome directly infer the **regulatory grammar** perturbed by non-coding variants. This approach enables a systematic and mechanistically grounded functional annotation of PGC3 risk loci, complementing existing proximity- and correlation-based methods.

## 2. Results

### 2.1. Identification and Classification of High-Confidence Drivers
To avoid the limitations of "Lead SNP" analysis and selection bias ("Winner's Curse"), we prioritized **10,832 variants** ($P < 5 \times 10^{-8}$) clustered into 237 loci using **Posterior Probabilities** derived from probabilistic fine-mapping.

To interpret the top hits biologically, we classified high-scoring genes into three distinct functional categories (Table 1).
*Note: Classification was based on brain expression specificity, functional annotation, and intolerance metrics (pLI), rather than statistical score alone.*

**Table 1: Functional Classification of Top Regulatory Drivers**
| Category | Criteria | Example Genes | Z-Score | Interpretation |
| :--- | :--- | :--- | :--- | :--- |
| **Regulatory Machinery** | Histones, TFs, Chromatin Modifiers | *H2AC20*, *ILF2*, *SETD1A* | 113.3 | Upstream modulators of gene expression networks. |
| **Neuronal Effectors** | Synaptic, Ion Channels, Transporters | *ATP2A2*, *SYNGAP1*, *ITPR3* | 16.5 | Downstream mediators of physiological function. |
| **Pleiotropic/Housekeeping** | Ubiquitous expression, High constitutive activity | *SERPINC1*, *RPS27* | 120.7 | Likely targets of high-density regulatory hubs; may represent broad cellular stress. |

This classification reveals that while "Regulatory Machinery" genes (like *H2AC20*) bear the highest statistical burden, the physiological specificity is provided by the "Neuronal Effectors" (*ATP2A2*).

### 2.2. Robustness of the Transcriptional Signal
The identification of *H2AC20* and *SERPINC1* as top hits raises the question of whether the Transcriptional Dysregulation signal is driven solely by these extreme outliers.
*   **Sensitivity Analysis:** To test robustness, we removed the top 1% highest-burden genes (including *H2AC20* and *SERPINC1*) and repeated the Unbiased GSEA.
*   **Result:** While the extreme statistical hyper-significance ($P < 10^{-30}$) attenuated, as expected given the removal of primary drivers, core biological modules including **Synapse Organization** ($P = 0.007$) and **Heterochromatin Organization** ($P = 0.01$) remained significant.
*   **Conclusion:** Notably, synaptic and chromatin-related terms remained enriched despite the removal of the most extreme regulatory hubs, indicating that convergence is not driven by a single locus but is a widely distributed polygenic property.

### 2.3. The Primary Point of Polygenic Convergence
In the unbiased GSEA against 5,000 biological processes, the strongest signal was **Negative Regulation of Transcription** ($P < 10^{-30}$). This suggests that the **primary point of polygenic convergence** in SCZ genetics is the **chromatin machinery** that governs gene expression. Histone genes should be interpreted as markers of this regulatory convergence rather than individual causal effectors.

### 2.4. Downstream Consequences: Synapse and Calcium
This upstream transcriptional dysregulation is accompanied by specific downstream failures:
*   **Vesicle-Mediated Transport** ($P = 8.4 \times 10^{-24}$)
*   **Intracellular Calcium Homeostasis** ($P = 3.5 \times 10^{-6}$): Validating the disruption of internal calcium handling (*ATP2A2*) as a specific mechanism.

### 2.5. Dominant Cellular Context
To resolve developmental timing, we applied **Footprint-Aware Normalization**.
*   **Result:** The bias towards progenitors vanished after normalization ($P>0.05$).
*   **Conclusion:** Risk is enriched in **Mature Excitatory Neurons** ($1.11x, P=0.007$) and **Glial Lineages** (Microglia: 1.32x, OPC: 1.39x), identifying a broader neuro-glial regulatory risk.

## 3. Discussion
We present a refined molecular architecture for Schizophrenia. This framework reconciles prior reports of synaptic, calcium, and epigenetic involvement by placing them within a single regulatory cascade:
1.  **Mechanism:** Polygenic risk disrupts **Enhancer Syntax** regulating chromatin modifiers (*H2AC20*). Here, "enhancer syntax" refers to the sequence-level organization of transcription factor binding motifs and spacing patterns learned implicitly by AlphaGenome, rather than explicit motif enumeration.
2.  **Propagation:** This transcriptional instability creates a failure of **Synaptic Maintenance** and **Calcium Buffering** (*ATP2A2*).
3.  **Context:** The pathology involves both **neuronal** and **glial** regulatory networks, broadening the scope beyond the previous "neuron-only" consensus.

We note that histone genes themselves (e.g., *H2AC20*) are unlikely to be dosage-sensitive causal drivers; rather, their recurrent identification reflects convergence on chromatin regulatory states that are broadly encoded by these loci. While neuronal regulatory burden remains dominant (affecting more loci), enrichment in microglia and oligodendrocyte progenitors suggests secondary modulation of risk via neuroimmune and myelination-related processes rather than primary causal drivers.

While AlphaGenome is a deep neural network, all downstream inferences in this study are based on calibrated, variant-level predictions aggregated via probabilistic fine-mapping rather than raw model scores. Future work will be required to resolve the precise directionality of regulatory effects at individual loci.

---

## 4. Methods
**Variant Prioritization:** PGC3 GWAS variants ($P < 5 \times 10^{-8}$) weighted by Approximate Bayes Factor Posterior Probabilities.
**Scoring:** AlphaGenome sequence-to-activity prediction. Scores represent the maximum absolute predicted log-fold change (LFC) across brain-relevant regulatory tracks.
**Robustness:** Sensitivity analysis performed by excluding top 1 percentile of gene scores.
**Enrichment:** Unbiased GSEA (Mann-Whitney U) against GO Biological Process 2023.
**Cell Type Normalization:** Binomial test adjusted for Total Accessible Base Pairs (bp) per cell type.

---

## 6. Figure Legends

**Figure 1: Deep Learning Decodes the Regulatory Syntax of Schizophrenia Risk.**
**(a)** Study design. 10,832 genome-wide significant variants ($P < 5 \times 10^{-8}$) from PGC3 were retained and clustered into 237 loci.
**(b)** Schematic of the **Posterior-Weighted Aggregation** method. Unlike standard "Lead SNP" approaches, AlphaGenome scores are weighted by the variant's approximate posterior probability ($PP$), ensuring that gene scores reflect the cumulative burden of likely causal variants while downweighting LD artifacts.
**(c)** Empirical Calibration. Gene scores are standardized against a whole-genome null distribution ($N=47,808$), yielding Z-scores that represent the deviation from expected background regulatory load.

**Figure 2: Polygenic Risk Converges on Transcriptional Machinery.**
**(a)** Regulatory Manhattan Plot. The gene-level regulatory burden (Z-score) across the genome. High-confidence drivers include *H2AC20* (Chr 6), *SERPINC1* (Chr 1), and *ILF2* (Chr 1).
**(b)** Unbiased GSEA. Rank-ordered enrichment of 5,000 GO Biological Processes. The top signal, **Negative Regulation of Transcription** ($P < 10^{-30}$), dominates the landscape.
**(c)** Sensitivity Analysis. Enrichment significance remains robust ($P < 0.01$) even after removing the top 1% of highest-scoring genes, confirming that this convergence is a broad polygenic property and not an artifact of outlier loci.

**Figure 3: A Hierarchy of Regulatory Dysfunction.**
**(a)** Functional Heatmap. Enriched pathways clustered by biological scaling. "Transcriptional Regulation" terms overlap with "Chromatin Assembly," while distinct clusters form for "Synaptic Transport" and "Calcium Homeostasis."
**(b)** The Calcium Module. Zoom-in on the **Intracellular Calcium Homeostasis** pathway, highlighting the specific disruption of *ATP2A2* (SERCA2) and *ITPR3* relative to voltage-gated channels.
**(c)** Regulatory-to-Effector Flow. Network schematic illustrating how upstream hits (Histones/TFs) logically precede downstream deficits (Vesicle Transport).

**Figure 4: Resolving the Cellular Context of Risk.**
**(a)** Artifactual Enrichment. Standard binomial testing (based on peak counts) suggests enrichment in iPSCs (stem cells).
**(b)** Normalization Logic. Correction for **Total Accessible Chromatin Volume** reveals that "iPSC enrichment" is a confounding effect of broader open chromatin.
**(c)** Definitive Specificity. After footprint-aware normalization, risk is enriched in **Mature Excitatory Neurons** ($1.11x$) and **Glial Cells** (Microglia/OPC > $1.3x$), overturning the simple "neuron-only" model.

**Figure 5: An Integrative Regulatory Framework.**
Conceptual synthesis of the findings. Polygenic non-coding variants constitute a "grammatical error" in the enhancer syntax of **neuronal chromatin regulators**. This leads to subtle but ubiquitous **transcriptional instability**, which disproportionately impacts systems with high metabolic and transport demands—specifically the **synaptic vesicle cycle** and **calcium buffering**. This reconciles the "Epigenetic" and "Synaptic" hypotheses into a single causal timeline.
