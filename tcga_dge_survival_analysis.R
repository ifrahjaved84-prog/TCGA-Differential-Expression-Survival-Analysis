# ==============================================================================
# PROJECT: TCGA RNA-Seq Differential Expression & Survival Analysis
# AUTHOR: Ifrah Javed Iqbal
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(survival)
  library(survminer)
  library(dplyr)
})

# 1. Generate Transcriptomic Counts & Clinical Survival Data
set.seed(42)
num_genes <- 1000
num_samples <- 60

count_matrix <- matrix(
  rnbinom(num_genes * num_samples, size = 2, mu = 200),
  nrow = num_genes, ncol = num_samples
)
rownames(count_matrix) <- paste0("Gene_", 1:num_genes)
colnames(count_matrix) <- paste0("Sample_", 1:num_samples)

# Elevate expression for candidate marker genes in tumor group
count_matrix[1:10, 1:30] <- count_matrix[1:10, 1:30] * 5

sample_info <- data.frame(
  sample_id = colnames(count_matrix),
  condition = factor(rep(c("Tumor", "Normal"), each = 30), levels = c("Normal", "Tumor")),
  days_to_death = sample(300:2000, num_samples, replace = TRUE),
  vital_status = sample(c(0, 1), num_samples, replace = TRUE, prob = c(0.4, 0.6))
)
rownames(sample_info) <- sample_info$sample_id

# 2. Differential Expression Analysis (DESeq2)
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_info, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))

res_df <- as.data.frame(res) %>%
  filter(!is.na(padj)) %>%
  mutate(Significance = case_when(
    log2FoldChange > 1 & padj < 0.05 ~ "Up-regulated",
    log2FoldChange < -1 & padj < 0.05 ~ "Down-regulated",
    TRUE ~ "Not Significant"
  ))

# 3. Generate Volcano Plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Up-regulated" = "firebrick", "Down-regulated" = "navy", "Not Significant" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Volcano Plot: Tumor vs Normal Expression", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value")

# 4. Kaplan-Meier Survival Modeling
vsd <- vst(dds, blind = FALSE)
top_gene <- head(rownames(res_df[order(res_df$padj), ]), 1)
expr_vals <- assay(vsd)[top_gene, sample_info$condition == "Tumor"]
tumor_samples <- sample_info[sample_info$condition == "Tumor", ]
tumor_samples$group <- ifelse(expr_vals > median(expr_vals), "High Expression", "Low Expression")

fit <- survfit(Surv(days_to_death, vital_status) ~ group, data = tumor_samples)

ggsurvplot(fit, data = tumor_samples, pval = TRUE, palette = c("firebrick", "dodgerblue4"),
           title = paste("Kaplan-Meier Curve:", top_gene), xlab = "Time (Days)")