# 🧬 Reference-Based RNA-Seq Data Analysis

> Splice-aware alignment, read quantification, differential expression analysis, and pathway enrichment for *Drosophila melanogaster* — run end-to-end on the **Galaxy platform**.

![Genome](https://img.shields.io/badge/Reference%20Genome-dm6-blue) ![Tool](https://img.shields.io/badge/Aligner-RNA%20STAR-purple) ![DE](https://img.shields.io/badge/DE%20Tool-DESeq2-red) ![Platform](https://img.shields.io/badge/Platform-Galaxy-brightgreen) ![Status](https://img.shields.io/badge/Status-Complete-brightgreen)

---

## 📋 Table of Contents

- [Project Overview](#-project-overview)
- [Dataset Overview](#-dataset-overview)
- [Pipeline Summary](#️-pipeline-summary)
- [01 — Alignment with STAR](#-01_alignment_star)
- [02 — Read Quantification with featureCounts](#-02_featurecounts)
- [03 — Differential Expression with DESeq2](#-03_deseq2_analysis)
- [04 — KEGG Pathway Visualization](#-04_kegg_pathways)
- [05 — GO Enrichment Analysis](#-05_go_enrichment)
- [Repository Structure](#️-repository-structure)
- [Tools & References](#️-tools--references)
- [Academic Context](#-academic-context)

---

## 🔬 Project Overview

This repository contains all outputs from a **reference-based RNA-Seq analysis pipeline** run on the [Galaxy platform](https://usegalaxy.org/), following the [RNA-Bio tutorial](https://rnabio.org/course/).

The dataset is a subset of the **Pasilla dataset** (Brooks et al., 2011), comparing *Drosophila melanogaster* samples with and without RNAi knockdown of the *pasilla* splicing factor gene. The pipeline covers the complete workflow — from raw read alignment through to pathway-level biological interpretation.

### This Repository Includes

| Asset | Description |
|---|---|
| 🗺️ STAR Alignment Logs | Per-sample mapping statistics from RNA STAR |
| 📊 featureCounts Summaries | Read-to-gene assignment tables |
| 📈 DESeq2 Outputs | Dispersion plots, MA plot, PCA, heatmaps |
| 🧭 KEGG Pathway Maps | Fold-change overlaid KEGG diagrams |
| 🔎 GO Enrichment Results | Top over-represented biological categories |

---

## 🦠 Dataset Overview

**Organism:** *Drosophila melanogaster*
**Reference Genome:** dm6 (Ensembl)
**Experimental Design:** Treated (*pasilla* RNAi knockdown) vs. Untreated

| Sample Accession | Condition | Library Type |
|---|---|---|
| GSM461176 | Untreated | Single-end |
| GSM461177 | Untreated | Paired-end |
| GSM461178 | Untreated | Paired-end |
| GSM461179 | Treated | Single-end |
| GSM461180 | Treated | Paired-end |
| GSM461181 | Treated | Paired-end |
| GSM461182 | Untreated | Single-end |

> Only **paired-end** samples (GSM461177 and GSM461180) are represented in the alignment logs and featureCounts summaries in this repository.

---

## ⚙️ Pipeline Summary

```
1. Raw Reads (FASTQ)
        ↓
2. Splice-Aware Alignment    →  RNA STAR        →  BAM files
        ↓
3. Read Quantification       →  featureCounts   →  Count matrices
        ↓
4. Differential Expression   →  DESeq2          →  DEG list + plots
        ↓
5. Pathway Visualization     →  Pathview        →  KEGG diagrams
        ↓
6. Functional Enrichment     →  goseq           →  GO term enrichment
```

---

## 📁 01_Alignment_STAR

**Tool:** RNA STAR *(Spliced Transcripts Alignment to a Reference)*

RNA STAR performs splice-aware alignment of raw reads to the *D. melanogaster* reference genome. It identifies exon-exon junctions using annotated splice sites (sjdb), making it well-suited for eukaryotic transcriptomics where multi-exon splicing is common.

### Files

| File | Description |
|---|---|
| `GSM461180_treat_paired_STAR_mapping_log.txt` | Alignment log — treated paired-end sample |
| `GSM461177_untreat_paired_STAR_mapping_log.txt` | Alignment log — untreated paired-end sample |

### Key Mapping Statistics

| Metric | GSM461180 (Treated) | GSM461177 (Untreated) |
|---|---|---|
| Input Reads | 1,042,655 | 1,116,426 |
| **Uniquely Mapped %** | **83.14%** | **78.98%** |
| Multi-mapped % | 5.50% | 4.54% |
| Unmapped — Too Short % | 5.74% | 11.11% |
| Mismatch Rate | 0.77% | 1.73% |
| Total Splices | 95,056 | 94,314 |
| Annotated Splices (sjdb) | 94,197 | 92,597 |

> Both samples show **high unique mapping rates (>78%)**, confirming good data quality. The majority of splices are canonical GT/AG, consistent with a well-annotated reference genome.

---

## 📁 02_FeatureCounts

**Tool:** featureCounts *(Subread package)*

featureCounts quantifies the number of reads overlapping each annotated genomic feature (gene/exon) using the aligned BAM files and reference GTF. Reads mapping to multiple loci or ambiguous features are excluded.

### Files

| File | Description |
|---|---|
| `GSM461180_treat_paired_featureCounts_summary.tabular` | Read assignment summary — treated paired-end sample |
| `GSM461177_untreat_paired_featureCounts_summary.tabular` | Read assignment summary — untreated paired-end sample |

### Assignment Summary

| Status | GSM461180 (Treated) | GSM461177 (Untreated) |
|---|---|---|
| ✅ Assigned | 843,496 | 825,225 |
| ❌ Unmapped | 183,960 | 118,428 |
| ❌ Multi-mapping | 273,323 | 324,121 |
| ❌ No Features | 14,778 | 20,088 |
| ❌ Ambiguity | 23,508 | 21,521 |

> A large proportion of reads are successfully assigned to annotated features in both samples. Multi-mapping reads are excluded to prevent ambiguous gene-level counts.

---

## 📁 03_DESeq2_Analysis

**Tool:** DESeq2 *(via Galaxy DESeq2 wrapper)*

DESeq2 models RNA-Seq count data using a **negative binomial distribution**, normalizes counts using size factors, estimates per-gene dispersion, and applies the Wald test to identify statistically significant differentially expressed genes (DEGs) between treated and untreated conditions.

### Files

| File | Description |
|---|---|
| `DESeq2_full_report.pdf` | Full multi-page report: PCA, sample distances, dispersion, p-value histogram, MA plot |
| `sample_clustering_heatmap.pdf` | Hierarchical clustering heatmap of all samples |
| `sample_to_sample_distances.pdf` | Euclidean distance matrix heatmap |
| `dispersion_estimates.png` | Per-gene dispersion estimates with fitted trend |
| `MA_plot_treated_vs_untreated.png` | Log2 fold change vs. mean expression for all genes |

### Output Descriptions

---

**🔷 PCA Plot** *(DESeq2_full_report.pdf — Page 1)*

Principal Component Analysis across all 7 samples. **PC1 (48% variance)** separates treated from untreated; **PC2 (33% variance)** separates paired-end (PE) from single-end (SE) libraries. This confirms that treatment effect is the dominant biological signal in the dataset.

---

**🔷 Sample-to-Sample Distance Heatmap** *(DESeq2_full_report.pdf — Page 2)*

A symmetric Euclidean distance matrix across all samples. Samples cluster correctly by condition and library type, confirming reproducible experimental grouping with no obvious outliers.

---

**🔷 Dispersion Estimates**

Per-gene dispersion estimates used to model count variance across samples. Black dots = per-gene maximum likelihood estimates; red line = fitted mean-dispersion trend; blue dots = final empirical Bayes-shrunken estimates used for testing. Dispersion decreasing with increasing mean expression is the expected and healthy pattern in RNA-Seq data.

![Dispersion Estimates](03_DESeq2_Analysis/dispersion_estimates.png)

---

**🔷 MA Plot**

Log2 fold change (Y-axis) plotted against mean normalized counts (X-axis, log10 scale). **Blue dots** = statistically significant DEGs (adjusted p-value < 0.05). A large number of significant genes are distributed across the full expression range, with fold changes reaching up to ±4 log2. Low-count genes (left side) show higher scatter, as expected.

![MA Plot — Treated vs Untreated](03_DESeq2_Analysis/MA_plot_treated_vs_untreated.png)

---

**🔷 p-value Histogram** *(DESeq2_full_report.pdf — Page 4)*

A sharp enrichment of p-values near 0 combined with a roughly uniform distribution across the remainder of the range confirms the presence of many true positives. This is the hallmark of a well-powered differential expression experiment.

---

## 📁 04_KEGG_Pathways

**Tool:** Pathview *(R package, via Galaxy)*

Pathview overlays DESeq2 log2 fold change values directly onto KEGG pathway diagrams for *D. melanogaster*. **Red boxes** = upregulated in treated; **green boxes** = downregulated in treated. This step contextualizes DEGs within established biological systems and metabolic networks.

### Files

| File | Pathway ID | Description |
|---|---|---|
| `dme00010_Glycolysis_Gluconeogenesis.png` | dme00010 | Glycolysis / Gluconeogenesis — metabolic enzyme expression changes |
| `dme03040_Spliceosome.png` | dme03040 | Spliceosome — expression changes in splicing machinery components |

### Pathway Interpretations

---

**🔷 Glycolysis / Gluconeogenesis (dme00010)**

Several glycolytic enzymes are differentially expressed upon *pasilla* knockdown, with notable upregulation in reactions involving glucose-6-phosphate and pyruvate metabolism. This suggests a secondary metabolic effect of *pasilla* loss on cellular energy homeostasis.

![Glycolysis / Gluconeogenesis Pathway](04_KEGG_Pathways/dme00010_Glycolysis_Gluconeogenesis.png)

---

**🔷 Spliceosome (dme03040)**

Since *pasilla* is a conserved RNA-binding splicing regulator, spliceosomal components are among the most biologically relevant targets. The pathway map highlights differential expression among U1/U2 snRNP-associated proteins (shown in red/green), directly consistent with *pasilla*'s known role in alternative splicing control.

![Spliceosome Pathway](04_KEGG_Pathways/dme03040_Spliceosome.png)

---

## 📁 05_GO_Enrichment

**Tool:** goseq *(R package, via Galaxy)*

Gene Ontology (GO) enrichment analysis identifies biological processes (BP), molecular functions (MF), and cellular components (CC) that are over-represented among the DEGs. goseq corrects for **gene-length bias** inherent to RNA-Seq using the Wallenius approximation, producing more reliable enrichment results than a standard hypergeometric test.

### Files

| File | Description |
|---|---|
| `GO_enrichment_top_categories.pdf` | Bubble plot — top over-represented GO categories across CC, BP, and MF ontologies |

### Top Enriched GO Categories

Bubble size = number of DEGs in that category; colour intensity = adjusted p-value (darker blue = more significant).

| GO Term | Ontology | % DE in Category | Interpretation |
|---|---|---|---|
| Glycogen biosynthetic process | BP | ~80% | Strongly altered carbohydrate storage metabolism |
| Glucan biosynthetic process | BP | ~80% | Overlaps with glycogen biosynthesis |
| Establishment of blood-brain barrier | BP | ~50% | Barrier and junction function affected |
| Septate junction assembly | BP | ~37% | Cell adhesion and paracellular sealing |
| Glutathione transferase activity | MF | ~23% | Detoxification and oxidative stress response |
| Transferase activity (alkyl) | MF | ~22% | Altered enzyme class activity |
| Small molecule metabolic process | BP | ~12% | Broad metabolic reprogramming |
| Oxidoreductase activity | MF | ~12% | Redox enzyme dysregulation |
| Response to stress | BP | ~12% | General cellular stress response |
| Extracellular region | CC | ~12% | Secreted / extracellular protein changes |

> The most significantly enriched terms span **carbohydrate metabolism** (glycogen/glucan biosynthesis) and **cellular organization** (septate junction, blood-brain barrier), revealing that *pasilla* knockdown has broad functional consequences well beyond alternative splicing alone.

---

## 🗂️ Repository Structure

```
RNA-Seq-Analysis/
│
├── README.md
│
├── 01_Alignment_STAR/
│   ├── GSM461180_treat_paired_STAR_mapping_log.txt
│   └── GSM461177_untreat_paired_STAR_mapping_log.txt
│
├── 02_FeatureCounts/
│   ├── GSM461180_treat_paired_featureCounts_summary.tabular
│   └── GSM461177_untreat_paired_featureCounts_summary.tabular
│
├── 03_DESeq2_Analysis/
│   ├── DESeq2_full_report.pdf
│   ├── sample_clustering_heatmap.pdf
│   ├── sample_to_sample_distances.pdf
│   ├── dispersion_estimates.png
│   └── MA_plot_treated_vs_untreated.png
│
├── 04_KEGG_Pathways/
│   ├── dme00010_Glycolysis_Gluconeogenesis.png
│   └── dme03040_Spliceosome.png
│
└── 05_GO_Enrichment/
    └── GO_enrichment_top_categories.pdf
```

### Folder Descriptions

**`01_Alignment_STAR/`** — RNA STAR alignment logs per sample. Report mapping rate, splice junction counts, mismatch rates, and multi-mapping statistics used to assess data quality before downstream analysis.

**`02_FeatureCounts/`** — Tab-delimited featureCounts summary tables showing how many reads were assigned, unmapped, multi-mapped, or ambiguously assigned for each sample.

**`03_DESeq2_Analysis/`** — Full DESeq2 differential expression outputs: QC plots (PCA, sample distances), model diagnostics (dispersion estimates), and result visualizations (MA plot, p-value histogram). Includes a combined PDF report.

**`04_KEGG_Pathways/`** — Pathview-generated KEGG pathway diagrams with DESeq2 log2 fold changes overlaid as a red-green colour gradient on enzyme/gene nodes.

**`05_GO_Enrichment/`** — goseq GO enrichment bubble plot showing the top over-represented biological categories among statistically significant DEGs, corrected for gene-length bias.

---

## 🛠️ Tools & References

| Step | Tool | Citation |
|---|---|---|
| Read Alignment | RNA STAR | Dobin et al., *Bioinformatics* (2013) |
| Read Quantification | featureCounts | Liao et al., *Bioinformatics* (2014) |
| Differential Expression | DESeq2 | Love et al., *Genome Biology* (2014) |
| Pathway Visualization | Pathview | Luo & Bhatt, *Bioinformatics* (2013) |
| GO Enrichment | goseq | Young et al., *Genome Biology* (2010) |
| Pipeline Platform | Galaxy | Afgan et al., *Nucleic Acids Res.* (2018) |

---

## 🎓 Academic Context

This project models a real-world RNA-Seq bioinformatics workflow and reflects core competencies in:

- **Splice-Aware Read Alignment** — handling exon-exon junctions in eukaryotic transcriptomes
- **Count-Based Quantification** — converting alignments to gene-level expression matrices
- **Statistical Differential Expression** — negative binomial modelling and multiple testing correction
- **Pathway-Level Interpretation** — connecting gene lists to known biological systems
- **Functional Enrichment** — identifying over-represented GO terms with length-bias correction

The resulting outputs are suitable for downstream workflows including integrative multi-omics studies, alternative splicing analysis, and biomarker discovery pipelines.

---

*Data sourced from [NCBI GEO — GSE18508](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508). Pipeline executed on [Galaxy](https://usegalaxy.org/) following the [RNA-Bio course](https://rnabio.org/course/).*