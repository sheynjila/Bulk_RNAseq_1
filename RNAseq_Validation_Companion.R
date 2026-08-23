# RNA-seq count/metadata validation companion script

validate_rnaseq_inputs <- function(counts, metadata, expected_samples = NULL) {
  if (!is.matrix(counts) && !is.data.frame(counts)) stop("counts must be a matrix or data frame.")
  if (!is.data.frame(metadata)) stop("metadata must be a data frame.")
  if (is.null(colnames(counts)) || is.null(rownames(metadata))) stop("Missing sample identifiers.")

  colnames(counts) <- trimws(colnames(counts))
  rownames(metadata) <- trimws(rownames(metadata))

  if (anyNA(colnames(counts)) || any(colnames(counts) == "")) stop("Empty count sample IDs.")
  if (anyNA(rownames(metadata)) || any(rownames(metadata) == "")) stop("Empty metadata sample IDs.")
  if (anyDuplicated(colnames(counts))) stop("Duplicate count sample IDs.")
  if (anyDuplicated(rownames(metadata))) stop("Duplicate metadata sample IDs.")

  count_only <- setdiff(colnames(counts), rownames(metadata))
  metadata_only <- setdiff(rownames(metadata), colnames(counts))
  if (length(count_only) || length(metadata_only)) {
    stop(paste0("Sample mismatch. Counts only: ", paste(count_only, collapse=", "),
                "; metadata only: ", paste(metadata_only, collapse=", ")))
  }

  if (!is.null(expected_samples) && !setequal(colnames(counts), trimws(expected_samples))) {
    stop("Observed samples do not match the expected sample list.")
  }

  metadata <- metadata[colnames(counts), , drop = FALSE]
  stopifnot(identical(colnames(counts), rownames(metadata)))
  if (!all(vapply(counts, is.numeric, logical(1)))) stop("Count columns must be numeric.")
  if (anyNA(counts)) stop("Counts contain missing values.")
  if (any(as.matrix(counts) < 0)) stop("Counts contain negative values.")

  message("Validation passed: ", ncol(counts), " samples are aligned.")
  list(counts=counts, metadata=metadata)
}

# Example import and repair workflow (edit paths before running):
# cnt <- read.csv("data/counts.csv", row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
# met <- read.csv("data/metadata.csv", check.names=FALSE, stringsAsFactors=FALSE)
# if ("X" %in% colnames(cnt)) cnt <- cnt[, colnames(cnt) != "X", drop=FALSE]
# rownames(met) <- trimws(as.character(met[[1]]))
# met <- met[, -1, drop=FALSE]
# validated <- validate_rnaseq_inputs(cnt, met)
# cnt <- validated$counts
# met <- validated$metadata
