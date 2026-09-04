# ==============================================================================
# ALIGN THE COUNT MATRIX AND SAMPLE METADATA
# Run this script from the project root after organizing the project files.
# ==============================================================================

data_directory <- "data"
counts_file <- file.path(data_directory, "counts.csv")
metadata_file <- file.path(data_directory, "metadata2.csv")

if (!file.exists(counts_file) || !file.exists(metadata_file)) {
  stop("counts.csv and metadata2.csv must be inside the data directory.")
}

# The supplied CSV files have a blank ID-column header stored as a trailing comma.
# Reading the first line separately keeps the first sample (SRR1039508) from being
# dropped or assigned the wrong counts.
sample_ids <- scan(
  counts_file,
  what = character(),
  sep = ",",
  nlines = 1,
  quiet = TRUE
)
sample_ids <- sample_ids[!is.na(sample_ids) & nzchar(sample_ids)]

# Read gene IDs as row names and skip the sample-name line already read above.
cnt <- read.csv(
  counts_file,
  header = FALSE,
  skip = 1,
  row.names = 1,
  check.names = FALSE
)

if (length(sample_ids) != ncol(cnt)) {
  stop("The number of sample names does not match the number of count columns.")
}
colnames(cnt) <- sample_ids

# The first metadata column contains sample IDs; the other two are variables.
met <- read.csv(
  metadata_file,
  header = FALSE,
  skip = 1,
  row.names = 1,
  check.names = FALSE
)

if (ncol(met) != 2) {
  stop("metadata2.csv must contain sample ID, cell line, and treatment columns.")
}
colnames(met) <- c("cellline", "dexamethasone")

# Stop early if IDs are duplicated or a count sample has no metadata row.
if (anyDuplicated(colnames(cnt))) {
  stop("Duplicate sample IDs were found in counts.csv.")
}
if (anyDuplicated(rownames(met))) {
  stop("Duplicate sample IDs were found in metadata2.csv.")
}

missing_metadata <- setdiff(colnames(cnt), rownames(met))
if (length(missing_metadata) > 0) {
  stop(
    "These count samples are missing from metadata2.csv: ",
    paste(missing_metadata, collapse = ", ")
  )
}

# Reorder metadata so its rows exactly match the count-matrix columns.
met <- met[colnames(cnt), , drop = FALSE]
stopifnot(identical(colnames(cnt), rownames(met)))

# Save files that can be loaded directly by LOAD_CLEANED_ALIGNED_DATA.R.
write.csv(
  cnt,
  file.path(data_directory, "counts_aligned.csv"),
  row.names = TRUE
)
write.csv(
  met,
  file.path(data_directory, "metadata_aligned.csv"),
  row.names = TRUE
)

cat("SUCCESS: aligned", ncol(cnt), "samples and", nrow(cnt), "genes.\n")
