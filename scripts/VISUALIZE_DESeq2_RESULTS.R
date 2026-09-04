# ==============================================================================
# VISUALIZE THE MAIN DESEQ2 RESULTS
# Run LOAD_CLEANED_ALIGNED_DATA.R before this script.
# ==============================================================================

library(DESeq2)
library(ggplot2)
library(dplyr)

if (!exists("dds") || !exists("res")) {
  stop("Run LOAD_CLEANED_ALIGNED_DATA.R before this script.")
}

padj_cutoff <- 0.05
lfc_cutoff <- 1

# Variance-stabilized counts are suitable for sample-level visualizations.
vsd <- vst(dds, blind = FALSE)

# ==============================================================================
# 1. PCA PLOT: DO SAMPLES CLUSTER BY TREATMENT OR CELL LINE?
# ==============================================================================

pca_data <- plotPCA(
  vsd,
  intgroup = c("dexamethasone", "cellline"),
  returnData = TRUE
)
percent_variance <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = dexamethasone, shape = cellline)
) +
  geom_point(size = 4, alpha = 0.8) +
  labs(
    title = "PCA: Sample Clustering",
    x = paste0("PC1: ", percent_variance[1], "% variance"),
    y = paste0("PC2: ", percent_variance[2], "% variance"),
    color = "Treatment",
    shape = "Cell Line"
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

print(pca_plot)

# ==============================================================================
# 2. VOLCANO PLOT: WHICH GENES CHANGE MOST STRONGLY?
# ==============================================================================

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

res_df <- res_df %>%
  mutate(
    significance = case_when(
      !is.na(padj) & padj < padj_cutoff & log2FoldChange > lfc_cutoff ~ "Up-regulated",
      !is.na(padj) & padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Down-regulated",
      TRUE ~ "Not significant"
    ),
    minus_log10_padj = -log10(pmax(padj, .Machine$double.xmin))
  )

volcano_plot <- res_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  ggplot(aes(x = log2FoldChange, y = minus_log10_padj, color = significance)) +
  geom_point(alpha = 0.6, size = 1.8) +
  scale_color_manual(
    values = c(
      "Up-regulated" = "#E41A1C",
      "Down-regulated" = "#377EB8",
      "Not significant" = "grey60"
    )
  ) +
  geom_vline(
    xintercept = c(-lfc_cutoff, lfc_cutoff),
    linetype = "dashed",
    color = "black",
    alpha = 0.5
  ) +
  geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype = "dashed",
    color = "black",
    alpha = 0.5
  ) +
  labs(
    title = "Volcano Plot: Dexamethasone vs Untreated",
    x = "Log2 fold change",
    y = "-Log10 adjusted p-value",
    color = NULL
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

print(volcano_plot)

# ==============================================================================
# 3. TOP-GENE BOXPLOT: HOW DOES THE STRONGEST RESULT LOOK PER SAMPLE?
# ==============================================================================

significant_genes <- res_df %>%
  filter(!is.na(padj), padj < padj_cutoff) %>%
  arrange(padj)

if (nrow(significant_genes) == 0) {
  message("No genes passed padj < ", padj_cutoff, "; the top-gene plot was skipped.")
} else {
  top_gene <- significant_genes$gene[1]

  gene_counts <- plotCounts(
    dds,
    gene = top_gene,
    intgroup = c("dexamethasone", "cellline"),
    returnData = TRUE
  )

  top_gene_plot <- ggplot(
    gene_counts,
    aes(x = dexamethasone, y = count + 1, fill = dexamethasone)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(aes(shape = cellline), width = 0.15, size = 3) +
    scale_y_log10() +
    labs(
      title = paste("Expression of Top Gene:", top_gene),
      x = "Treatment",
      y = "Normalized count + 1 (log10 scale)",
      fill = "Treatment",
      shape = "Cell Line"
    ) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  print(top_gene_plot)
}
