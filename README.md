# 🧬 Reference-Based RNA-Seq Data Analysis

> A complete walkthrough of an RNA-Seq differential gene expression pipeline on *Drosophila melanogaster* — from raw FASTQ reads to functional pathway enrichment, run end-to-end on the **Galaxy platform**.

![Genome](https://img.shields.io/badge/Reference%20Genome-dm6-blue) ![Aligner](https://img.shields.io/badge/Aligner-RNA%20STAR-purple) ![DE](https://img.shields.io/badge/DE%20Tool-DESeq2-red) ![Platform](https://img.shields.io/badge/Platform-Galaxy-brightgreen) ![Status](https://img.shields.io/badge/Status-Complete-brightgreen)

---

## 📋 Table of Contents

- [Project Overview](#-project-overview)
- [Biological Background](#-biological-background)
- [Dataset & Experimental Design](#-dataset--experimental-design)
- [Pipeline Summary](#️-pipeline-summary)
- [Step 1 — Quality Control](#-step-1--quality-control)
- [Step 2 — Spliced Mapping to Reference Genome](#-step-2--spliced-mapping-to-reference-genome)
- [Step 3 — Library Strandness Estimation](#-step-3--library-strandness-estimation)
- [Step 4 — Read Counting per Gene](#-step-4--read-counting-per-gene)
- [Step 5 — Differential Gene Expression Analysis](#-step-5--differential-gene-expression-analysis)
- [Step 6 — Functional Enrichment Analysis](#-step-6--functional-enrichment-analysis)
- [Key Results](#-key-results)
- [Repository Structure](#️-repository-structure)
- [File Formats Used](#-file-formats-used)
- [Tools & References](#️-tools--references)

---

## 🔬 Project Overview

This repository documents the **systematic identification of differentially expressed (DE) genes** between *Pasilla (PS) gene-depleted* and normal *Drosophila melanogaster* cells using RNA-Seq data. Starting from raw sequencing reads, the analysis performs quality filtering, splice-aware genome mapping, read quantification, statistical DE testing, and functional enrichment to understand the biological impact of PS depletion.

### This Repository Includes

| Asset | Description |
|---|---|
| 🗺️ STAR Alignment Logs | Per-sample mapping statistics |
| 📊 featureCounts Summaries | Read-to-gene assignment tables |
| 📈 DESeq2 Outputs | Dispersion plots, MA plot, PCA, heatmaps |
| 🧭 KEGG Pathway Maps | Fold-change overlaid KEGG diagrams |
| 🔎 GO Enrichment Results | Top over-represented biological categories |

**Organism:** *Drosophila melanogaster*
**Reference Genome:** dm6
**Original Study:** Brooks et al. 2011 — *Conservation of an RNA regulatory map between Drosophila and mammals*
**Tutorial:** [Galaxy Training Network RNA-Seq tutorial](https://gxy.io/GTN:T00295)

---

## 🦠 Biological Background

The **Pasilla (PS)** gene encodes a splicing regulator in *Drosophila*, homologous to the mammalian **Nova-1** and **Nova-2** proteins. In the original study, PS was depleted by RNA interference (RNAi) to identify which genes and molecular pathways it regulates.

| Condition | Description |
|---|---|
| **Treated** | Pasilla gene knocked down by RNAi |
| **Untreated** | Normal Pasilla expression |

Comparing RNA-Seq reads between these two conditions reveals which genes change expression when PS is absent — validating the experiment if the *pasilla* gene itself appears downregulated in treated samples.

---

## 🧫 Dataset & Experimental Design

| Sample ID | Condition | Library Type |
|---|---|---|
| GSM461176 | Untreated | Single-end |
| GSM461177 | Untreated | Paired-end |
| GSM461178 | Untreated | Paired-end |
| GSM461179 | Treated (PS depleted) | Single-end |
| GSM461180 | Treated (PS depleted) | Paired-end |
| GSM461181 | Treated (PS depleted) | Paired-end |
| GSM461182 | Untreated | Single-end |

**4 untreated replicates, 3 treated replicates.** Two treated and two untreated samples are paired-end; the rest are single-end. This mix is accounted for as a **second factor** in the DESeq2 statistical model to remove technical variation from library preparation.

---

## ⚙️ Pipeline Summary

```
Raw FASTQ reads
      │
      ▼
 Quality Control          ←  Falco, MultiQC, Cutadapt
      │
      ▼
 Spliced Mapping          ←  STAR (dm6 genome + GTF annotation)
      │
      ▼
 Mapping QC               ←  MultiQC (STAR logs), IGV / JBrowse2
      │
      ▼
 Strandness Estimation    ←  Infer Experiment (RSeQC)
      │
      ▼
 Read Counting            ←  featureCounts
      │
      ▼
 Differential Expression  ←  DESeq2 (multi-factor model)
      │
      ▼
 Annotation & Filtering   ←  DESeq2 annotator, Filter tool
      │
      ▼
 Visualization            ←  heatmap2 (normalized counts + Z-scores)
      │
      ▼
 Functional Enrichment    ←  goseq (GO + KEGG), Pathview
```

---

## 🔷 Step 1 — Quality Control

Raw reads contain sequencing errors and adapter contamination. This step removes low-quality bases and short/poor-quality reads **before mapping** to prevent misalignments that would corrupt downstream results.

### 1.1 Per-Sample Quality Check

**Tool:** Falco *(FastQC-compatible quality report generator)*

**What to check in the report:**
- Per-base sequence quality (should be mostly green / Q > 28)
- Adapter content (flags if adapters are present)
- Sequence length distribution
- Read length (36 bp for this dataset)

### 1.2 Aggregated QC Report

**Tool:** MultiQC — combines all Falco reports into a single HTML summary

**What to look for:**
- Quality drops toward the 3′ end of reads (common in Illumina data)
- Overrepresented sequences or adapter contamination
- Consistent quality across all samples

### 1.3 Trimming and Filtering

**Tool:** Cutadapt

| Parameter | Value | Reason |
|---|---|---|
| Quality cutoff (R1) | 20 | Remove bases with Phred score < 20 |
| Minimum read length | 20 bp | Discard reads too short to map reliably |
| Input type | Paired-end collection | Processes both mates together |

> **Note:** Paired-end samples are trimmed together to ensure consistency — if one read in a pair is discarded, its partner is too.

---

## 🔷 Step 2 — Spliced Mapping to Reference Genome

### Why Spliced Mapping?

In eukaryotes, mature mRNA has introns removed (spliced out). A standard aligner would fail to map reads that span two exons because the sequence doesn't exist as-is in the genomic DNA. **Splice-aware aligners** identify exon-exon junction reads and correctly map them across intron gaps.

```
Genome:    ████ exon1 ████----intron----████ exon2 ████
Read:                    ████████████████
                         (spans junction → needs spliced mapping)
```

### 2.1 Genome Mapping

**Tool:** STAR *(Spliced Transcripts Alignment to a Reference)*

| Parameter | Value | Reason |
|---|---|---|
| Reference genome | dm6 Full | Drosophila reference |
| GTF annotation | BDGP6.32.109 | Defines known splice sites |
| Junction overhang length | 36 | Read length − 1 |
| Per gene read counts | GeneCounts | Generates counts alongside mapping |
| Coverage output | Bedgraph (strand 1 & 2) | Used for strandness estimation |

**Output files:**

| File | Description |
|---|---|
| `mapped.bam` | Aligned reads in binary format |
| `log` | Mapping statistics (used for QC) |
| `reads per gene` | Raw gene-level counts (unstranded, forward, reverse columns) |
| `Coverage strand 1 & 2` | Per-base coverage bedgraphs |

### 2.2 Mapping Quality Check

**Tool:** MultiQC on STAR logs

| Metric | Expected | Action if not met |
|---|---|---|
| % Uniquely mapped | ~80% | < 70% → investigate contamination |
| % Multi-mapped | < 10% | High → check genome version match |
| % Unmapped | < 20% | High → check trimming or wrong genome |

#### Alignment Results for This Dataset

| Metric | GSM461180 (Treated PE) | GSM461177 (Untreated PE) |
|---|---|---|
| Input Reads | 1,042,655 | 1,116,426 |
| **Uniquely Mapped %** | **83.14%** | **78.98%** |
| Multi-mapped % | 5.50% | 4.54% |
| Unmapped — Too Short % | 5.74% | 11.11% |
| Mismatch Rate | 0.77% | 1.73% |
| Total Splices | 95,056 | 94,314 |
| Annotated Splices (sjdb) | 94,197 | 92,597 |

> Both samples achieve **>78% unique mapping rates** — confirming good data quality and correct genome alignment. The majority of splices are canonical GT/AG.

### 2.3 Visual Inspection of Alignments

**Tools:** IGV, Sashimi Plot, JBrowse2 — inspect BAM alignments at `chr4:540,000–560,000`

**What to observe:**
- **Grey coverage peaks** — read depth across the region; higher = more expression
- **Connecting arcs between reads** — spliced reads spanning introns; arc = skipped intron
- **Sashimi plots** — arcs with numbers represent splice junctions; number = supporting reads

---

## 🔷 Step 3 — Library Strandness Estimation

### Why Strandness Matters

Some library preparation protocols preserve information about which strand the original RNA came from. This affects how reads are counted per gene — particularly for reads in regions where two genes on opposite strands overlap. Using the wrong strandness setting inflates or deflates counts for many genes.

| Library Type | Description |
|---|---|
| Unstranded | Reads map to genes on both strands equally (~0.5 / ~0.5) |
| Stranded forward | Reads map predominantly to the same strand as the gene |
| Stranded reverse | Reads map predominantly to the opposite strand of the gene |

**Tool:** Infer Experiment (RSeQC) — samples 200,000 reads and reports the fraction explained by each strand orientation.

> **Result for this dataset:** Both fractions ≈ 0.5 → **Unstranded library**

This strandness setting is carried forward into featureCounts.

---

## 🔷 Step 4 — Read Counting per Gene

### Purpose

Quantify how many reads (fragments) map to each annotated gene, producing the raw count matrix used for differential expression testing.

**Tool:** featureCounts

| Parameter | Value |
|---|---|
| Strandness | Unstranded (determined in Step 3) |
| Feature type | exon |
| Gene identifier | gene_id |
| Count paired reads as | 1 fragment |
| Minimum mapping quality | 10 |

**Output files:**

| File | Description |
|---|---|
| `Counts` | Tab-delimited: Gene ID → read count per sample |
| `Feature lengths` | Gene ID → length in bp (required for goseq) |
| `Summary` | Assignment statistics (% assigned vs unassigned) |

#### featureCounts Assignment Summary

| Status | GSM461180 (Treated) | GSM461177 (Untreated) |
|---|---|---|
| ✅ Assigned | 843,496 | 825,225 |
| ❌ Unmapped | 183,960 | 118,428 |
| ❌ Multi-mapping | 273,323 | 324,121 |
| ❌ No Features | 14,778 | 20,088 |
| ❌ Ambiguity | 23,508 | 21,521 |

> An assignment rate below 50% is a warning sign — causes include genome/annotation version mismatch, incorrect strandness, or poor mapping quality.

### Why Raw Counts Can't Be Compared Directly

| Bias | Cause | Effect |
|---|---|---|
| Sequencing depth | Different samples have different total reads | Deeper samples have higher counts for all genes |
| Gene length | Longer genes produce more reads | Long genes appear more expressed |
| Library composition | One condition may express unique high-abundance genes | Remaining genes appear falsely downregulated |

DESeq2 corrects for all of these biases during normalization.

---

## 🔷 Step 5 — Differential Gene Expression Analysis

**Tool:** DESeq2 *(negative binomial model with multi-factor design)*

DESeq2 normalizes count data using size factors, estimates per-gene dispersion using empirical Bayes shrinkage, and applies the Wald test to identify statistically significant DEGs.

### Experimental Design — Two Factors

| Factor | Levels |
|---|---|
| Treatment *(primary)* | treated → GSM461179, 180, 181 / untreated → GSM461176, 177, 178, 182 |
| Sequencing type *(covariate)* | PE → GSM461177, 178, 180, 181 / SE → GSM461176, 179, 182 |

> Including sequencing type as a second factor removes technical variation from library preparation method, producing a cleaner treatment effect estimate.

### DESeq2 Result Table Columns

| Column | Description |
|---|---|
| GeneID | Flybase gene ID (e.g., FBgn0003360) |
| Base mean | Mean normalized count across all 7 samples |
| log₂(FC) | Log₂ fold change: treated vs untreated (positive = up in treated) |
| StdErr | Standard error of the log₂FC estimate |
| Wald-Stat | Test statistic |
| P-value | Raw p-value |
| **P-adj** | **Benjamini-Hochberg adjusted p-value (FDR-corrected) — use this** |

### Diagnostic Plot Interpretations

---

**🔹 Dispersion Estimates**

Black dots = per-gene maximum likelihood dispersion estimates. Red line = fitted mean-dispersion trend. Blue dots = final empirical Bayes-shrunken estimates used for hypothesis testing. Dispersion decreasing with increasing mean expression is the expected and healthy pattern — it confirms the model fit is reliable.

![Dispersion Estimates](03_DESeq2_Analysis/dispersion_estimates.png)

---

**🔹 MA Plot**

X-axis: mean normalized counts (log10 scale). Y-axis: log₂ fold change. **Blue dots** = statistically significant DEGs (adjusted p-value < 0.05). Significant genes are spread across all expression levels, with fold changes up to ±4 log₂. Low-count genes (left) show higher scatter — this is expected and biologically uninformative.

![MA Plot — Treated vs Untreated](03_DESeq2_Analysis/MA_plot_treated_vs_untreated.png)

---

**🔹 PCA Plot** *(DESeq2_full_report.pdf — Page 1)*

PC1 (48% variance) separates **treated vs untreated** — confirming the treatment has a strong and dominant effect. PC2 (33% variance) separates **paired-end vs single-end** libraries — confirming the second factor was necessary and effective.

---

**🔹 Sample-to-Sample Distance Heatmap** *(DESeq2_full_report.pdf — Page 2)*

Symmetric Euclidean distance matrix across all 7 samples. Treated samples cluster together; untreated cluster together — confirming reproducible experimental grouping with no obvious outliers.

---

**🔹 p-value Histogram** *(DESeq2_full_report.pdf — Page 4)*

A sharp enrichment near p = 0, combined with a roughly uniform distribution for the rest, confirms the presence of many true positives. This is the hallmark of a well-powered differential expression experiment with a genuine biological signal.

---

### Filtering DE Genes

Two sequential filters applied to annotated DESeq2 results:

| Filter | Condition | Rationale |
|---|---|---|
| Statistical significance | Adjusted p-value < 0.05 | Controls false discovery rate at 5% |
| Biological relevance | \|log₂FC\| > 1 (i.e., > 2-fold change) | Removes statistically significant but negligibly small changes |

**Result:** ~113 genes pass both filters.

### Heatmap of Normalized Counts (Top DE Genes)

113 DE genes × 7 samples, clustered by expression profile. Treated samples cluster together; untreated cluster together. Expression patterns split cleanly along the treatment axis, visually confirming the DE results.

Log₂(value + 1) transformation applied to stabilize variance for visualization.

---

### Z-score Heatmap

Same 113 genes, but Z-score normalized per gene across samples. Z-score removes the effect of absolute expression level and highlights **relative patterns** — making it easier to see which samples are high or low for each gene compared to the average.

```
Z-score formula:  z = (x_ij − mean_i) / sd_i
```

---

## 🔷 Step 6 — Functional Enrichment Analysis

### 6.1 Gene Ontology (GO) Enrichment

**Tool:** goseq — specifically corrects for **gene-length bias** inherent to RNA-Seq data (longer genes have more reads and are more likely to be called DE, which inflates GO term enrichment without correction).

**Input 1:** Boolean DE gene list (True/False per gene based on adjusted p-value < 0.05)
**Input 2:** Gene length file from featureCounts

| Parameter | Value |
|---|---|
| Genome | Fruit fly (dm6) |
| Gene ID format | Ensembl Gene ID |
| Categories | GO: Cellular Component, Biological Process, Molecular Function |

**Output table key columns:**

| Column | Description |
|---|---|
| category | GO term ID |
| numDEInCat | Number of DE genes in this category |
| numInCat | Total genes annotated to this category |
| term | GO term description |
| ontology | BP / MF / CC |
| **p.adjust.over_represented** | **BH-adjusted p-value — use this for significance** |

### 6.2 KEGG Pathway Analysis

**Tool:** goseq (same inputs, KEGG category selected)

**Notable pathways identified:**

| KEGG ID | Pathway | Result |
|---|---|---|
| dme00010 | Glycolysis / Gluconeogenesis | **Over-represented** (adj. p < 0.05) |
| dme03040 | Spliceosome | **Under-represented** (biologically expected: PS is a splicing regulator) |

### 6.3 KEGG Pathway Visualization

**Tool:** Pathview — overlays log₂FC values from DESeq2 onto KEGG pathway diagrams

**Colour legend:**
- 🟥 **Red box** = upregulated in treated samples
- 🟩 **Green box** = downregulated in treated samples
- ⬜ **Grey box** = gene present but not significantly DE
- ◻️ **White box** = gene not detected in dataset

---

**🔹 Glycolysis / Gluconeogenesis (dme00010)**

Several glycolytic enzymes are differentially expressed upon *pasilla* knockdown. Notable upregulation in reactions involving glucose-6-phosphate and pyruvate metabolism suggests a secondary effect on cellular energy homeostasis.

![Glycolysis / Gluconeogenesis Pathway](04_KEGG_Pathways/dme00010_Glycolysis_Gluconeogenesis.png)

---

**🔹 Spliceosome (dme03040)**

Since *pasilla* is a conserved RNA-binding splicing regulator, spliceosomal components are among the most biologically relevant targets. Differential expression among U1/U2 snRNP-associated proteins is directly consistent with *pasilla*'s known role in alternative splicing control.

![Spliceosome Pathway](04_KEGG_Pathways/dme03040_Spliceosome.png)

---

### Top Enriched GO Categories

Bubble size = number of DE genes in category; colour = adjusted p-value (darker = more significant). See `05_GO_Enrichment/GO_enrichment_top_categories.pdf` for the full plot.

| GO Term | Ontology | % DE in Category | Interpretation |
|---|---|---|---|
| Glycogen biosynthetic process | BP | ~80% | Strongly altered carbohydrate storage metabolism |
| Glucan biosynthetic process | BP | ~80% | Overlaps with glycogen biosynthesis |
| Establishment of blood-brain barrier | BP | ~50% | Barrier/junction function affected |
| Septate junction assembly | BP | ~37% | Cell adhesion and paracellular sealing |
| Glutathione transferase activity | MF | ~23% | Detoxification and oxidative stress response |
| Transferase activity (alkyl) | MF | ~22% | Altered enzyme class activity |
| Small molecule metabolic process | BP | ~12% | Broad metabolic reprogramming |
| Oxidoreductase activity | MF | ~12% | Redox enzyme dysregulation |
| Response to stress | BP | ~12% | General cellular stress response |
| Extracellular region | CC | ~12% | Secreted/extracellular protein changes |

> The most significant terms span **carbohydrate metabolism** and **cellular organization**, revealing that *pasilla* knockdown has broad consequences beyond alternative splicing alone.

---

## 📊 Key Results

| Analysis Stage | Result |
|---|---|
| Read mapping rate | ~80% uniquely mapped (both test samples) |
| Library strandness | Unstranded |
| Total genes tested | ~14,000 |
| Genes with adj p < 0.05 | Several hundred |
| Genes with adj p < 0.05 AND \|log₂FC\| > 1 | **~113 genes** |
| Pasilla gene (FBgn0261552) | Downregulated in treated ✅ *(validates experiment)* |
| Top over-represented GO terms | RNA binding, splicing, nucleic acid processing |
| Significant KEGG pathway | dme00010 — Glycolysis / Gluconeogenesis |

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

**`01_Alignment_STAR/`** — RNA STAR alignment logs reporting mapping rate, splice junction counts, mismatch rates, and multi-mapping statistics. Used to assess data quality before downstream quantification.

**`02_FeatureCounts/`** — Tab-delimited assignment summary tables showing how many reads were assigned, unmapped, multi-mapped, or ambiguously assigned per sample.

**`03_DESeq2_Analysis/`** — Full DESeq2 outputs: QC plots (PCA, sample distances heatmap), model diagnostics (dispersion estimates), and result visualizations (MA plot, p-value histogram). Combined PDF report included.

**`04_KEGG_Pathways/`** — Pathview-generated KEGG pathway diagrams with DESeq2 log₂FC values overlaid as a red-green colour gradient on enzyme/gene nodes.

**`05_GO_Enrichment/`** — goseq GO enrichment bubble plot showing top over-represented biological categories among significant DEGs, corrected for gene-length bias.

---

## 📄 File Formats Used

| Format | Extension | Description |
|---|---|---|
| FASTQ | `.fastqsanger` | Raw sequencing reads + quality scores |
| BAM | `.bam` | Binary file of read-to-genome alignments |
| GTF | `.gtf.gz` | Gene annotation (exon positions, gene IDs) |
| BED12 | `.bed` | 12-column annotation format (for Infer Experiment) |
| Counts | `.counts` / `.tabular` | Tab-delimited: Gene ID → integer read count |
| Bedgraph | `.bedgraph` | Per-base read coverage across the genome |

---

## 🛠️ Tools & References

| Step | Tool | Citation |
|---|---|---|
| QC Reports | Falco | — |
| Aggregate QC | MultiQC | Ewels et al., *Bioinformatics* (2016) |
| Adapter Trimming | Cutadapt | Marcel, *EMBnet.journal* (2011) |
| Read Alignment | RNA STAR | Dobin et al., *Bioinformatics* (2013) |
| Strandness | Infer Experiment (RSeQC) | Wang et al., *Bioinformatics* (2012) |
| Read Quantification | featureCounts | Liao et al., *Bioinformatics* (2014) |
| Differential Expression | DESeq2 | Love et al., *Genome Biology* (2014) |
| Heatmaps | heatmap2 | — |
| GO & KEGG Enrichment | goseq | Young et al., *Genome Biology* (2010) |
| Pathway Visualization | Pathview | Luo & Brouwer, *Bioinformatics* (2013) |
| Pipeline Platform | Galaxy | Afgan et al., *Nucleic Acids Res.* (2018) |

**Dataset:** Brooks AN et al. (2011). *Conservation of an RNA regulatory map between Drosophila and mammals.* Genome Research 21:193–202. [GEO: GSE18508](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18508)

---

*Pipeline completed following the [Galaxy Training Network RNA-Seq tutorial](https://gxy.io/GTN:T00295) (CC-BY 4.0).*