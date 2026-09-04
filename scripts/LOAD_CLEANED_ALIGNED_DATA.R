# ==============================================================================
# LOAD THE ALIGNED DATA AND RUN DESEQ2
# Run Aligning_meta_counts.R once before running this script.
# ==============================================================================

library(DESeq2)

data_directory <- "data"

# Load gene counts and sample metadata. Column 1 contains row names in both files.
cnt <- read.csv(
  file.path(data_directory, "counts_aligned.csv"),
  row.names = 1,
  check.names = FALSE
)
met <- read.csv(
  file.path(data_directory, "metadata_aligned.csv"),
  row.names = 1,
  check.names = FALSE
)

required_columns <- c("cellline", "dexamethasone")
if (!all(required_columns %in% colnames(met))) {
  stop(
    "metadata_aligned.csv must contain: ",
    paste(required_columns, collapse = ", "),
    ". Run Aligning_meta_counts.R again."
  )
}

# DESeq2 requires count columns and metadata rows in exactly the same order.
if (!identical(colnames(cnt), rownames(met))) {
  stop("Sample names or sample order do not match between the two aligned files.")
}

count_matrix <- as.matrix(cnt)
if (!is.numeric(count_matrix) || anyNA(count_matrix)) {
  stop("The count matrix must contain only numeric values and no missing values.")
}
if (any(count_matrix < 0) || any(count_matrix != round(count_matrix))) {
  stop("DESeq2 counts must be non-negative whole numbers.")
}

# Set treatment and cell line as categorical variables.
met$cellline <- factor(met$cellline)
met$dexamethasone <- factor(met$dexamethasone)

required_treatments <- c("untreated", "treated")
if (!all(required_treatments %in% levels(met$dexamethasone))) {
  stop("The treatment column must contain both 'untreated' and 'treated'.")
}
met$dexamethasone <- relevel(met$dexamethasone, ref = "untreated")

# The design adjusts the treatment comparison for differences among cell lines.
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = met,
  design = ~ cellline + dexamethasone
)

# Remove genes with fewer than 10 total reads across all samples.
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Fit the model and explicitly compare treated with untreated samples.
dds <- DESeq(dds)
res <- results(
  dds,
  contrast = c("dexamethasone", "treated", "untreated"),
  alpha = 0.05
)

cat("SUCCESS: DESeq2 analysis completed for", nrow(dds), "genes.\n")
summary(res)
head(res[order(res$padj), ])
