# TCGA RNA-Seq Differential Expression & Survival Analysis Pipeline

An end-to-end computational pipeline in R designed to identify differentially expressed genes (DEGs) in tumor vs. normal tissue and evaluate candidate markers using Kaplan-Meier overall survival modeling. 

This project bridges clinical pharmacy principles with translational bioinformatics, serving as a framework for target discovery, biomarker validation, and precision oncology research.

---

## Technical Stack & Dependencies

* **Language:** R (v4.6+)
* **Core Libraries:**
  * `DESeq2` — Differential gene expression (DGE) statistical modeling
  * `survival` & `survminer` — Kaplan-Meier survival curves and log-rank statistics
  * `ggplot2` — Publication-ready data visualizations (Volcano plots)
  * `dplyr` & `pheatmap` — Data manipulation and heatmapping

---

## Pipeline Overview

1. **Synthetic RNA-Seq Dataset Generation:**
   * Simulated transcriptomic count matrix representing 1,000 genes across 60 samples (30 Tumor, 30 Normal) using a negative binomial distribution (`rnbinom`).
   * Matched clinical metadata including follow-up time (`days_to_death`) and patient `vital_status`.

2. **Differential Expression Analysis (`DESeq2`):**
   * Filtered low-count reads ($\sum \text{counts} \ge 10$).
   * Calculated $\log_2$ fold changes and Benjamini-Hochberg adjusted $p$-values (`padj`).
   * Categorized genes into **Up-regulated** ($\log_2\text{FC} > 1, \text{padj} < 0.05$), **Down-regulated** ($\log_2\text{FC} < -1, \text{padj} < 0.05$), or **Not Significant**.

3. **Prognostic Survival Modeling:**
   * Applied variance stabilizing transformations (`vst`) to normalize transcriptomic counts.
   * Stratified tumor cohorts into **High Expression** and **Low Expression** groups relative to median baseline expression.
   * Modeled overall survival probability over time using the Kaplan-Meier estimator and log-rank tests.

---

## Key Visualizations

### 1. Differential Expression (Volcano Plot)
Identifies statistically significant transcriptomic alterations between tumor and normal patient cohorts.

![Volcano Plot](volcano_plot.png)

### 2. Kaplan-Meier Survival Analysis
Evaluates whether a candidate overexpressed gene serves as a statistically significant prognostic driver for overall survival ($p = 0.49$).

![Kaplan-Meier Curve](km_curve.png)

---

## Biological & Clinical Takeaway

While `Gene_1` showed significant differential expression between tumor and normal samples during transcriptomic screening, subsequent Kaplan-Meier survival modeling yielded $p = 0.49$. 

**Translational Insight:** Overexpression in tumor tissue alone is insufficient to declare a gene a functional driver or therapeutic drug target. Cross-referencing differential expression with patient survival metrics is essential to distinguish true driver mutations from passenger biomarkers in rational target discovery.

---

## Author

**Ifrah **  
*Doctor of Pharmacy (PharmD) | Computational Oncology & Precision Medicine Research*
