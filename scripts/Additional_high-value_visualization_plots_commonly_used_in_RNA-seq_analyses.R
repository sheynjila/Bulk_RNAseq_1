# ==============================================================================
# ADDITIONAL RNA-SEQ PLOTS
# Run LOAD_CLEANED_ALIGNED_DATA.R before this script.
# VISUALIZE_DESeq2_RESULTS.R does not need to be run first.
# ==============================================================================

# Install packages only once, from the R console, if they are missing:
# install.packages(c("ggplot2", "dplyr", "tidyr", "pheatmap", "RColorBrewer"))
# BiocManager::install(c(
#   "DESeq2", "apeglm", "clusterProfiler", "enrichplot",
#   "org.Hs.eg.db", "pathview"
# ))

library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

if (!exists("dds") || !exists("res")) {
  stop("Run LOAD_CLEANED_ALIGNED_DATA.R before this script.")
}

# Change an option to TRUE only when you want to run that optional section.
run_lfc_shrinkage <- FALSE
run_enrichment <- FALSE
run_pathview <- FALSE

padj_cutoff <- 0.05

# Create the two objects used throughout this script.
vsd <- vst(dds, blind = FALSE)
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

significant_results <- res_df %>%
  filter(!is.na(padj), padj < padj_cutoff) %>%
  arrange(padj)

# ==============================================================================
# 1. SAMPLE-TO-SAMPLE DISTANCE HEATMAP
# Large distances or unexpected clusters can reveal outliers or batch effects.
# ==============================================================================

sample_distances <- dist(t(assay(vsd)))
sample_distance_matrix <- as.matrix(sample_distances)
distance_colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(
  sample_distance_matrix,
  clustering_distance_rows = sample_distances,
  clustering_distance_cols = sample_distances,
  color = distance_colors,
  main = "Sample-to-Sample Distance Matrix"
)

# ==============================================================================
# 2. TOP-30 SIGNIFICANT-GENE HEATMAP
# Row Z-scores show each gene's relative expression pattern across samples.
# ==============================================================================

top_genes <- head(significant_results$gene, 30)

if (length(top_genes) == 0) {
  message("No genes passed padj < ", padj_cutoff, "; the gene heatmap was skipped.")
} else {
  heatmap_matrix <- assay(vsd)[top_genes, , drop = FALSE]

  # Constant rows cannot be converted to Z-scores, so remove them if present.
  row_sd <- apply(heatmap_matrix, 1, sd)
  heatmap_matrix <- heatmap_matrix[
    is.finite(row_sd) & row_sd > 0,
    ,
    drop = FALSE
  ]

  if (nrow(heatmap_matrix) == 0) {
    message("The selected genes had no variation; the gene heatmap was skipped.")
  } else {
    scaled_matrix <- t(scale(t(heatmap_matrix)))
    sample_annotation <- as.data.frame(
      colData(vsd)[, c("cellline", "dexamethasone")]
    )

    pheatmap(
      scaled_matrix,
      annotation_col = sample_annotation,
      cluster_rows = nrow(scaled_matrix) > 1,
      show_colnames = TRUE,
      show_rownames = TRUE,
      main = "Top Significant Genes (Row Z-score)"
    )
  }
}

# ==============================================================================
# 3. MA PLOT
# This shows fold change across the range of average gene expression.
# ==============================================================================

ma_plot <- res_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange), baseMean > 0) %>%
  mutate(significant = padj < padj_cutoff) %>%
  ggplot(aes(x = baseMean, y = log2FoldChange, color = significant)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_x_log10() +
  scale_color_manual(
    values = c("FALSE" = "grey60", "TRUE" = "#E41A1C"),
    labels = c("Not significant", "Significant")
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  labs(
    title = "MA Plot: Fold Change vs Mean Expression",
    x = "Mean normalized count (log10 scale)",
    y = "Log2 fold change",
    color = NULL
  ) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

print(ma_plot)

# ==============================================================================
# 4. RAW P-VALUE DISTRIBUTION
# A peak near zero suggests real signal; the remaining values are often flatter.
# ==============================================================================

pvalue_plot <- res_df %>%
  filter(!is.na(pvalue)) %>%
  ggplot(aes(x = pvalue)) +
  geom_histogram(bins = 50, fill = "#377EB8", color = "white") +
  labs(
    title = "Distribution of Raw P-values",
    x = "P-value",
    y = "Number of genes"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

print(pvalue_plot)

# ==============================================================================
# 5. DESEQ2 DISPERSION DIAGNOSTIC
# The fitted trend should follow the center of the gene-level estimates.
# ==============================================================================

plotDispEsts(dds, main = "DESeq2 Dispersion Estimates")

# ==============================================================================
# 6. TOP-GENE EXPRESSION PROFILES
# Lines connect the mean expression for untreated and treated samples per cell line.
# ==============================================================================

top_10_genes <- head(significant_results$gene, 10)

if (length(top_10_genes) == 0) {
  message("No genes passed padj < ", padj_cutoff, "; the profile plot was skipped.")
} else {
  normalized_counts <- counts(dds, normalized = TRUE)[
    top_10_genes,
    ,
    drop = FALSE
  ]

  counts_long <- as.data.frame(normalized_counts) %>%
    tibble::rownames_to_column("gene") %>%
    pivot_longer(
      cols = -gene,
      names_to = "sample",
      values_to = "normalized_count"
    )

  sample_data <- as.data.frame(colData(dds)) %>%
    tibble::rownames_to_column("sample")

  profile_data <- counts_long %>%
    left_join(sample_data, by = "sample") %>%
    group_by(gene, cellline, dexamethasone) %>%
    summarise(
      mean_log2_expression = mean(log2(normalized_count + 1)),
      .groups = "drop"
    )

  profile_plot <- ggplot(
    profile_data,
    aes(
      x = dexamethasone,
      y = mean_log2_expression,
      group = gene,
      color = gene
    )
  ) +
    geom_line(alpha = 0.7) +
    geom_point(size = 2) +
    facet_wrap(~ cellline) +
    labs(
      title = "Top-Gene Expression Profiles",
      x = "Treatment",
      y = "Mean log2(normalized count + 1)",
      color = "Gene"
    ) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  print(profile_plot)
}

# ==============================================================================
# OPTIONAL 1. LOG2 FOLD-CHANGE SHRINKAGE
# Set run_lfc_shrinkage <- TRUE above after installing apeglm.
# Shrinkage reduces noisy fold changes, especially for low-count genes.
# ==============================================================================

if (run_lfc_shrinkage) {
  library(apeglm)

  treatment_coefficient <- grep(
    "^dexamethasone_",
    resultsNames(dds),
    value = TRUE
  )

  if (length(treatment_coefficient) != 1) {
    stop(
      "Could not identify one treatment coefficient. Available names: ",
      paste(resultsNames(dds), collapse = ", ")
    )
  }

  res_lfc <- lfcShrink(
    dds,
    coef = treatment_coefficient,
    type = "apeglm"
  )

  par(mfrow = c(1, 2))
  plotMA(res, ylim = c(-4, 4), main = "Unshrunken Fold Change")
  plotMA(res_lfc, ylim = c(-4, 4), main = "Shrunken Fold Change")
  par(mfrow = c(1, 1))
}

# ==============================================================================
# OPTIONAL 2. GO ENRICHMENT AND GSEA
# Set run_enrichment <- TRUE above after installing the listed Bioconductor packages.
# org.Hs.eg.db is the annotation package for human genes.
# ==============================================================================

if (run_enrichment) {
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)

  # Remove optional Ensembl version suffixes such as ".12".
  significant_ensembl <- unique(sub("\\..*$", "", significant_results$gene))

  if (length(significant_ensembl) == 0) {
    message("No significant genes were available for GO enrichment.")
  } else {
    go_results <- enrichGO(
      gene = significant_ensembl,
      OrgDb = org.Hs.eg.db,
      keyType = "ENSEMBL",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05
    )

    if (nrow(as.data.frame(go_results)) == 0) {
      message("No GO Biological Process terms passed the enrichment cutoff.")
    } else {
      go_plot <- dotplot(go_results, showCategory = 15) +
        ggtitle("GO Enrichment: Biological Processes") +
        theme_bw()
      print(go_plot)
    }
  }

  # GSEA ranks all tested genes, so it does not require a DEG cutoff.
  gsea_gene_ids <- sub("\\..*$", "", res_df$gene)
  gene_list <- res_df$stat
  names(gene_list) <- gsea_gene_ids
  gene_list <- gene_list[is.finite(gene_list) & !duplicated(names(gene_list))]
  gene_list <- sort(gene_list, decreasing = TRUE)

  gsea_results <- gseGO(
    geneList = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "BP",
    pvalueCutoff = 0.05,
    verbose = FALSE
  )

  gsea_table <- as.data.frame(gsea_results)
  if (nrow(gsea_table) == 0) {
    message("No GO gene sets passed the GSEA cutoff.")
  } else {
    print(
      gseaplot2(
        gsea_results,
        geneSetID = gsea_table$ID[1],
        title = gsea_table$Description[1]
      )
    )
  }
}

# ==============================================================================
# OPTIONAL 3. KEGG PATHWAY OVERLAY
# Set run_pathview <- TRUE above after installing pathview, clusterProfiler,
# and org.Hs.eg.db. Change pathway_id to the human KEGG pathway you need.
# ==============================================================================

if (run_pathview) {
  library(pathview)
  library(clusterProfiler)
  library(org.Hs.eg.db)

  pathway_id <- "hsa04110" # Cell cycle

  fold_changes <- res_df %>%
    transmute(
      ENSEMBL = sub("\\..*$", "", gene),
      log2FoldChange = log2FoldChange
    ) %>%
    filter(!is.na(log2FoldChange))

  gene_ids <- bitr(
    unique(fold_changes$ENSEMBL),
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  mapped_table <- inner_join(gene_ids, fold_changes, by = "ENSEMBL")
  mapped_fold_changes <- tapply(
    mapped_table$log2FoldChange,
    mapped_table$ENTREZID,
    mean
  )

  pathview(
    gene.data = mapped_fold_changes,
    pathway.id = pathway_id,
    species = "hsa",
    kegg.dir = "plots",
    limit = list(gene = 2, cpd = 1)
  )
}

# UpSet plots are not included here because this project currently has one
# treatment contrast. They become useful after two or more contrasts are analyzed.
