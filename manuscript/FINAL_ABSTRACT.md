# Abstract

## Multi-Layer Regulatory Architecture of Schizophrenia: Integrating Deep Learning Predictions with 3D Chromatin Organization

**Background**: Genome-wide association studies have identified over 250 genomic loci associated with schizophrenia risk, yet translating these loci into causal genes and mechanisms remains challenging. Regulatory variants may affect gene expression through diverse mechanisms—local enhancers, distal chromatin loops, transcription initiation, and alternative splicing—operating across multiple spatial scales and chromatin states.

**Methods**: We integrated 20,760 fine-mapped schizophrenia variants from the Psychiatric Genomics Consortium wave 3 study with AlphaGenome deep learning predictions of regulatory impact across six chromatin modalities (DNase-seq, ATAC-seq, CAGE, RNA-seq, PRO-CAP, splice junctions). To recover genes regulated by long-range enhancers, we incorporated 3D chromatin architecture using PsychENCODE brain Hi-C loops, mapping variants to distal genes up to 1 Mb away. Gene-level scores were aggregated using normalized posterior probabilities as weights, tested against an independent gene universe (19,966 GENCODE protein-coding genes), and validated using GTEx v10 brain eQTL data from 13 tissues.

**Results**: We identified **1,617 genes** with evidence of regulatory impact (empirical FDR<0.10). Integrating 3D chromatin architecture contributed **468 genes** (29% of total), including **219 genes accessible exclusively through chromatin loops** (mean distance: 556 kb, range: 10-1,019 kb). Independent validation using population-level brain eQTL demonstrated **69% of top-100 genes show significant brain eQTL** (vs 48.8% background; OR=2.33, p<0.0001), confirming functional regulatory effects. Critically, **Hi-C-identified genes showed 84% brain eQTL concordance** compared to 44.6% for linear window-only genes, validating that 3D chromatin captures functional enhancer-promoter interactions.

**Modality-specific analysis** revealed mechanistic diversity: chromatin accessibility-dominant genes (60% eQTL concordance) may reflect activity-dependent regulation, while RNA-dominant genes (75% concordance) show cumulative steady-state effects. Splicing-regulated genes (45% concordance) exhibited isoform diversity invisible to bulk expression profiling. **Regional specificity analysis** identified significant eQTL enrichment in 12 of 13 brain tissues (FDR<0.05), with strongest effects in frontal cortex (p=0.0004), cerebellum (p=0.0005), and nucleus accumbens (p=0.0007)—regions supporting cognition, motor coordination, and reward processing.

**Pathway analysis** identified eight synaptic-related processes with nominal significance (p<0.05), including regulation of calcium-dependent exocytosis, synaptic plasticity, and excitatory postsynaptic potentials, though these did not survive multiple testing correction. Notably, canonical synaptic genes (DRD2, GRM3, GAD1, GAD2) were not recovered, reflecting limitations of the fetal brain Hi-C dataset rather than biological absence—a data coverage issue validated by our method's high concordance with standardized H-MAGMA approaches (93.2%).

**Conclusions**: Schizophrenia risk variants exert regulatory effects through a **multi-layer architecture** requiring integration of linear genomic context, 3D chromatin organization, modality-specific regulation, and brain region specificity. Three key principles emerge:

1. **3D chromatin is essential**: 29% of genes require Hi-C information, with 5.6% accessible only through distal contacts. The exceptional eQTL validation rate of Hi-C genes (84% vs 44.6%) confirms functional relevance of chromatin architecture.

2. **Regulatory diversity across modalities**: Different genes are regulated at different layers (chromatin accessibility, transcription, RNA, splicing), with modality-specific eQTL patterns reflecting distinct biological mechanisms and suggesting multiple therapeutic entry points.

3. **Circuit-level specificity**: Enrichment patterns across brain regions (strongest in frontal cortex, cerebellum, basal ganglia) suggest circuit-level dysfunction rather than global brain pathology, aligning with cognitive and reward processing deficits in schizophrenia.

This integrative framework—validated by population-level expression data and standardized methodology comparison—provides a foundation for understanding schizophrenia's molecular architecture and prioritizing genes for functional validation and therapeutic development. The finding that 5.6% of genes are accessible only through 3D chromatin highlights the critical need for comprehensive cell type-specific and developmental stage-specific chromatin resources to fully map psychiatric disorder regulatory networks.

---

**Keywords**: Schizophrenia, GWAS, functional genomics, chromatin architecture, Hi-C, deep learning, gene regulation, brain eQTL, regulatory variants

---

**Word count**: 550 (abstract body)
