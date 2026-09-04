# ==============================================================================
# PACKAGE MANAGEMENT FOR THE RNA-SEQ PROJECT
#
# Run this script once to create the manage_packages() function. It lists package
# status by default. To use another mode, run one of these commands in R:
#
# manage_packages("list")       # List packages and installation status
# manage_packages("manual")     # Print commands for manual installation
# manage_packages("automatic")  # Install only packages that are missing
# ==============================================================================

# Packages used by the cleaned analysis and visualization scripts.
project_packages <- data.frame(
  package = c(
    "DESeq2",
    "ggplot2",
    "dplyr",
    "tidyr",
    "pheatmap",
    "RColorBrewer",
    "apeglm",
    "clusterProfiler",
    "enrichplot",
    "org.Hs.eg.db",
    "pathview"
  ),
  source = c(
    "Bioconductor",
    rep("CRAN", 5),
    rep("Bioconductor", 5)
  ),
  category = c(
    rep("Core", 6),
    rep("Optional advanced analysis", 5)
  ),
  purpose = c(
    "Differential-expression analysis",
    "General plotting",
    "Data manipulation",
    "Reshaping data for plots",
    "Heatmaps",
    "Heatmap color palettes",
    "Log2 fold-change shrinkage",
    "GO enrichment and GSEA",
    "Enrichment plots",
    "Human gene annotations",
    "KEGG pathway diagrams"
  ),
  stringsAsFactors = FALSE
)

manage_packages <- function(action = "list") {
  action <- match.arg(tolower(action), c("list", "manual", "automatic"))

  # Show the current status before doing anything else.
  package_status <- project_packages
  package_status$installed <- vapply(
    package_status$package,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )

  cat("\nRNA-seq project packages\n")
  print(package_status, row.names = FALSE)

  if (action == "list") {
    cat("\nInstalled:", sum(package_status$installed),
        "of", nrow(package_status), "packages.\n")
    return(invisible(package_status))
  }

  cran_packages <- project_packages$package[
    project_packages$source == "CRAN"
  ]
  bioconductor_packages <- project_packages$package[
    project_packages$source == "Bioconductor"
  ]

  if (action == "manual") {
    quote_packages <- function(packages) {
      paste0('"', packages, '"', collapse = ", ")
    }

    cat(
      "\nCopy and run these commands in the R console:\n\n",
      "# Install CRAN packages\n",
      "install.packages(c(", quote_packages(cran_packages), "))\n\n",
      "# Install the Bioconductor package manager if needed\n",
      "if (!requireNamespace(\"BiocManager\", quietly = TRUE)) {\n",
      "  install.packages(\"BiocManager\")\n",
      "}\n\n",
      "# Install Bioconductor packages\n",
      "BiocManager::install(c(",
      quote_packages(bioconductor_packages),
      "), ask = FALSE, update = FALSE)\n",
      sep = ""
    )

    return(invisible(package_status))
  }

  # Automatic mode installs only packages that are currently missing.
  missing_cran <- cran_packages[
    !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  missing_bioconductor <- bioconductor_packages[
    !vapply(bioconductor_packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_cran) > 0) {
    message("Installing missing CRAN packages: ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran)
  } else {
    message("All CRAN packages are already installed.")
  }

  if (length(missing_bioconductor) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }

    message(
      "Installing missing Bioconductor packages: ",
      paste(missing_bioconductor, collapse = ", ")
    )
    BiocManager::install(
      missing_bioconductor,
      ask = FALSE,
      update = FALSE
    )
  } else {
    message("All Bioconductor packages are already installed.")
  }

  still_missing <- project_packages$package[
    !vapply(
      project_packages$package,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(still_missing) == 0) {
    cat("\nSUCCESS: All project packages are installed.\n")
  } else {
    warning(
      "These packages are still missing: ",
      paste(still_missing, collapse = ", ")
    )
  }

  invisible(still_missing)
}

# Safe default: report package status without installing anything.
manage_packages("list")
