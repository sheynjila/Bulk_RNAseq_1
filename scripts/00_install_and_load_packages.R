# ==============================================================================
# INSTALL AND LOAD PACKAGES FOR THIS RNA-SEQ ANALYSIS
#
# This file demonstrates two approaches:
#   1. Manual installation: each package is written out separately.
#   2. Automated installation: R finds and installs only missing packages.
#
# Sourcing this file defines the functions below but does not install anything.
# ==============================================================================

# ==============================================================================
# METHOD 1: MANUAL INSTALLATION AND LOADING
# ==============================================================================

install_packages_manually <- function() {
  # CRAN packages
  install.packages("ggplot2")
  install.packages("dplyr")
  install.packages("tidyr")
  install.packages("pheatmap")
  install.packages("RColorBrewer")

  # BiocManager is needed to install Bioconductor packages.
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  # Bioconductor packages
  BiocManager::install("DESeq2", ask = FALSE, update = FALSE)
  BiocManager::install("apeglm", ask = FALSE, update = FALSE)
  BiocManager::install("clusterProfiler", ask = FALSE, update = FALSE)
  BiocManager::install("enrichplot", ask = FALSE, update = FALSE)
  BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)
  BiocManager::install("pathview", ask = FALSE, update = FALSE)
}

load_packages_manually <- function() {
  # Core analysis and plotting packages
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)

  # Optional advanced-analysis packages
  library(apeglm)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(pathview)

  cat("All RNA-seq packages were loaded successfully.\n")
}

# To use the manual method, run these commands in the R console:
# install_packages_manually()
# load_packages_manually()

# ==============================================================================
# METHOD 2: AUTOMATED INSTALLATION AND LOADING
# ==============================================================================

cran_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "pheatmap",
  "RColorBrewer"
)

bioconductor_packages <- c(
  "DESeq2",
  "apeglm",
  "clusterProfiler",
  "enrichplot",
  "org.Hs.eg.db",
  "pathview"
)

all_analysis_packages <- c(cran_packages, bioconductor_packages)

install_packages_automatically <- function() {
  # Identify CRAN packages that are not installed.
  missing_cran <- cran_packages[
    !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_cran) > 0) {
    message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran)
  } else {
    message("All required CRAN packages are already installed.")
  }

  # Identify Bioconductor packages that are not installed.
  missing_bioconductor <- bioconductor_packages[
    !vapply(
      bioconductor_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing_bioconductor) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }

    message(
      "Installing Bioconductor packages: ",
      paste(missing_bioconductor, collapse = ", ")
    )
    BiocManager::install(
      missing_bioconductor,
      ask = FALSE,
      update = FALSE
    )
  } else {
    message("All required Bioconductor packages are already installed.")
  }
}

load_packages_automatically <- function() {
  missing_packages <- all_analysis_packages[
    !vapply(
      all_analysis_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing_packages) > 0) {
    stop(
      "These packages are not installed: ",
      paste(missing_packages, collapse = ", "),
      ". Run install_packages_automatically() first."
    )
  }

  invisible(
    lapply(
      all_analysis_packages,
      library,
      character.only = TRUE
    )
  )

  cat("All RNA-seq packages were loaded successfully.\n")
}

# To use the automated method, run these commands in the R console:
# install_packages_automatically()
# load_packages_automatically()

cat(
  "Package functions are ready. See METHOD 1 or METHOD 2 in this script.\n"
)
