# Sequence-Based Deep Learning Identifies ATP2A2-Centered Calcium Dysregulation in Schizophrenia Polygenic Risk

**Authors:** Yusuf Cicek, et al.
**Affiliation:** Istanbul University - Cerrahpasa

## Abstract
Schizophrenia (SCZ) GWAS discovery has outpaced functional understanding, identifying hundreds of non-coding loci whose mechanisms remain obscure. Here, we applied **AlphaGenome**, a deep learning framework, to decode the regulatory syntax of SCZ risk. Initially targeting 257 lead SNPs, we identified a putative convergence on intracellular calcium signaling (*ATP2A2*, *ITPR3*). To rigorously validate this finding against linkage disequilibrium (LD) structure and selection biases, we expanded our analysis to 10,832 candidate causal variants using a **Posterior-Weighted Credible Set** approach. This comprehensive, unbiased modeling revealed that SCZ risk is fundamentally driven by a massive dysregulation of **neuronal transcriptional control**, which causally propagates to specific deficits in **synaptic transport** and **intracellular calcium homeostasis**. By implementing chromatin-footprint normalization, we further resolved that **risk enrichment extends beyond neurons to include microglia and OPCs**, challenging the simple "neuron-only" consensus. These findings support a coherent regulatory framework linking non-coding risk to transcriptional and synaptic pathways.

...

## 2. Phase I: The Primary Screen (Lead SNPs)
...
*   **Discovery:** Gene Set Enrichment Analysis (GSEA) identified a significant enrichment of **Intracellular Calcium Signaling** ($P = 0.009$).
*   **Key Candidates:** The SERCA2 calcium pump (*ATP2A2*) and IP3 Receptor (*ITPR3*) emerged as top targets, **identifying disruption of internal calcium handling as a predicted mechanism**.

...

## 4. Discussion
We present a refined **regulatory landscape** for Schizophrenia. This framework reconciles prior reports of synaptic, calcium, and epigenetic involvement by placing them within a single regulatory cascade:
1.  **Mechanism:** ...
2.  **Propagation:** ...
3.  **Context:** The pathology involves both **neuronal** and **glial** regulatory networks, **consistent with emerging evidence for neuro-glial interactions in SCZ pathogenesis**.

We note that histone genes themselves (e.g., *H2AC20*) are unlikely to be dosage-sensitive causal drivers; rather, their recurrent identification reflects convergence on chromatin regulatory states that are broadly encoded by these loci. While neuronal regulatory burden remains dominant (affecting more loci), enrichment in microglia and oligodendrocyte progenitors suggests secondary modulation of risk via neuroimmune and myelination-related processes rather than primary causal drivers.

**Limitations:**
Our study has several limitations. First, as with any sequence-based predictor, AlphaGenome captures regulatory potential rather than realized transcriptional output, and thus complements rather than replaces experimental assays. Second, while we prioritized variants using fine-mapping, the precise causal nucleotide at each locus requires experimental confirmation. Third, our cell-type specificity analysis relies on publicly available ATAC-seq atlases, which may not capture all developmental stages. Finally, while we identify a convergence on Calcium signaling, the specific downstream physiological consequences await in vivo validation.

## 5. Methods Summary
*   **GWAS:** PGC3 Schizophrenia ($N=320k$).
*   **Model:** AlphaGenome (DeepMind API).
*   **Scoring:** AlphaGenome sequence-to-activity prediction. Scores represent the maximum absolute predicted log-fold change (LFC) across **53 brain-relevant regulatory tracks**.
...

---

## 1. Introduction
Genome-wide association studies (GWAS) have identified 287 genomic loci associated with schizophrenia (Trubetskoy et al., 2022). However, standard post-GWAS methods (H-MAGMA, proximity mapping) often rely on broad chromatin correlations that lack the resolution to distinguish causal regulatory variants from LD proxies.

We deployed **AlphaGenome**, a sequence-based deep learning model, to predict the cell-type-specific impact of variants on gene expression and chromatin topology. This approach enables a systematic and mechanistically grounded functional annotation of PGC3 risk loci, complementing existing proximity- and correlation-based methods.

## 2. Phase I: The Primary Screen (Lead SNPs)
Our initial study focused on the 257 lead SNPs defined by the PGC3 GWAS.
*   **Method:** We scored each lead variant for regulatory impact across 53 relevant biosamples.
*   **Discovery:** Gene Set Enrichment Analysis (GSEA) identified a significant enrichment of **Intracellular Calcium Signaling** ($P = 0.009$).
*   **Key Candidates:** The SERCA2 calcium pump (*ATP2A2*) and IP3 Receptor (*ITPR3*) emerged as top targets.
*   **Validation:** We observed strong overlap with official PGC3 prioritized genes (OR=5.03, $P < 10^{-10}$), providing initial confidence in the model.

## 3. Phase II: Methodological Refinement & Rigor
To elevate this finding to a definitive "causal" model, we addressed critical limitations inherent in "Lead SNP" analyses (LD blindness, selection bias) through a series of rigorous methodological upgrades.

### 3.1. Credible Set Prioritization (Addressing "LD Blindness")
Lead SNPs often tag lower-impact causal variants. We implemented a **Proxy Credible Set** strategy:
*   **Expansion:** We extracted all variants with $P < 5 \times 10^{-8}$ ($N=10,832$) clustered into 237 loci.
*   **Weighting:** We calculated "Proxy Posterior Probabilities" (PP) using Approximate Bayes Factors to weight variants based on their statistical likelihood of causality.

### 3.2. AlphaGenome Scoring & Aggregation (Addressing "Winner's Curse")
To avoid overestimating risk from noisy variants, we replaced simple "max scores" with **Posterior-Weighted Aggregation**:
$$GeneScore = \sum (AlphaGenome\_LFC \times Proxy\_PP)$$
This ensures that only variants with both *high regulatory impact* and *high causal probability* drive the gene score.

### 3.3. Unbiased Whole-Genome GSEA (Addressing "Universe Bias")
Our initial Calcium finding was hypothesis-driven. To test its robustness, we performed an **Unbiased Whole-Genome GSEA** against 5,000 Biological Processes ($N=47,808$ genes).
*   **Result 1 (The Core):** The strongest signal was **Negative Regulation of Transcription** ($P < 10^{-30}$), driven by histone modifiers (*H2AC20*) and transcriptional regulators (*ILF2*).
*   **Result 2 (The mechanism):** The **Intracellular Calcium Homeostasis** signal survived the unbiased competition ($P < 10^{-5}$), confirming it as a specific downstream mechanism.
*   **Result 3 ( The Output):** These disruptions converge on **Synaptic Vesicle Transport** ($P < 10^{-24}$).

### 3.4. Cellular Specificity (Addressing "Footprint Bias")
Previous studies implicated fetal progenitors (iPSCs). We suspected this was confounded by the large open chromatin volume in stem cells.
*   **Method:** We applied a **Footprint-Aware Binomial Test**, normalizing expected hits by Total Accessible Base Pairs.
*   **Result:** iPSC enrichment vanished ($P=0.07$).
*   **Conclusion:** Significant specific enrichment persists in **Mature Excitatory Neurons** (1.11x, $P < 0.01$) and **Glial Lineages** (Microglia/OPC ~ 1.35x), identifying a broad regulatory landscape.

---

## 4. Discussion
We present a refined **regulatory landscape** for Schizophrenia. This framework reconciles prior reports of synaptic, calcium, and epigenetic involvement by placing them within a single regulatory cascade:
1.  **Mechanism:** Polygenic risk disrupts **Enhancer Syntax** regulating chromatin modifiers (*H2AC20*). Here, "enhancer syntax" refers to the sequence-level organization of transcription factor binding motifs and spacing patterns learned implicitly by AlphaGenome, rather than explicit motif enumeration.
2.  **Propagation:** This transcriptional instability creates a failure of **Synaptic Maintenance** and **Calcium Buffering** (*ATP2A2*).
3.  **Context:** The pathology involves both **neuronal** and **glial** regulatory networks, **refining the previously neuron-dominant model**.

We note that histone genes themselves (e.g., *H2AC20*) are unlikely to be dosage-sensitive causal drivers; rather, their recurrent identification reflects convergence on chromatin regulatory states that are broadly encoded by these loci. While neuronal regulatory burden remains dominant (affecting more loci), enrichment in microglia and oligodendrocyte progenitors suggests secondary modulation of risk via neuroimmune and myelination-related processes rather than primary causal drivers.

As with any sequence-based predictor, AlphaGenome captures regulatory potential rather than realized transcriptional output, and thus complements rather than replaces experimental assays. While AlphaGenome is a deep neural network, all downstream inferences in this study are based on calibrated, variant-level predictions aggregated via probabilistic fine-mapping rather than raw model scores. Future work will be required to resolve the precise directionality of regulatory effects at individual loci.

---

## 5. Methods Summary
*   **GWAS:** PGC3 Schizophrenia ($N=320k$).
*   **Model:** AlphaGenome (DeepMind API).
*   **Scoring:** AlphaGenome sequence-to-activity prediction. Scores represent the maximum absolute predicted log-fold change (LFC) across brain-relevant regulatory tracks.
*   **Statistics:** P-Weighted Aggregation, Empirical Z-Scoring ($N=47k$), FDR-Corrected GSEA (Mann-Whitney U).
*   **Code Availability:** All scripts are available in the `scz_hypothesis_testing/` directory.
