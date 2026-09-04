## =============================================================================
## install_rnaseq_packages.R
##
## Installs the R/Bioconductor packages needed for a full RNA-seq workflow,
## from raw FASTQ reads through to differential expression and downstream
## enrichment analysis. Organism focus: human (hg38 / Ensembl).
##
## Workflow stage -> package group mapping:
##   1. Bootstrap                  -> BiocManager
##   2. Raw read QC & trimming     -> Rfastp, ShortRead, Biostrings
##   3. Alignment (R-native)       -> Rsubread, Rhisat2, Rbowtie2, QuasR
##   4. BAM handling & quant.      -> Rsamtools, GenomicAlignments, GenomicFeatures,
##                                     GenomicRanges, SummarizedExperiment,
##                                     tximport, tximeta
##   5. Human annotation & genome  -> org.Hs.eg.db, TxDb.Hsapiens.UCSC.hg38.knownGene,
##                                     EnsDb.Hsapiens.v86, AnnotationDbi, AnnotationHub,
##                                     ensembldb, biomaRt, BSgenome.Hsapiens.UCSC.hg38
##   6. Differential expression    -> DESeq2, apeglm, ashr, IHW, vsn
##   7. Downstream / enrichment    -> pheatmap, EnhancedVolcano, PCAtools,
##                                     clusterProfiler, enrichplot, DOSE, ReactomePA
##   8. Data wrangling / reporting -> tidyverse, rmarkdown, knitr
##
## Usage:
##   Rscript install_rnaseq_packages.R
##   (or source() it inside an R session)
## =============================================================================

## ---- 0. Bootstrap: BiocManager ---------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

message(sprintf("Using Bioconductor version: %s", as.character(BiocManager::version())))

## Helper: install any package (CRAN or Bioconductor) only if not already
## installed. BiocManager::install() transparently handles both sources.
install_if_missing <- function(pkgs) {
    for (pkg in pkgs) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            message(sprintf("  Installing %s ...", pkg))
            BiocManager::install(pkg, update = FALSE, ask = FALSE)
        } else {
            message(sprintf("  %s already installed - skipping.", pkg))
        }
    }
}

## ---- 1. Raw read QC & trimming (R-only, no external binaries needed) ------

message("\n[Section] Raw read QC & trimming ...")
qc_trimming_pkgs <- c(
    "Rfastp",     # bundles fastp binary: adapter trimming + QC report, all in R
    "ShortRead",  # FASTQ I/O and manual quality inspection
    "Biostrings"  # sequence handling used throughout the workflow
)
install_if_missing(qc_trimming_pkgs)

## ---- 2. Alignment (R-native aligners, multiple approaches) ----------------

message("\n[Section] Alignment ...")
alignment_pkgs <- c(
    "Rsubread", # built-in align()/subjunc() aligner + featureCounts()
    "Rhisat2",  # HISAT2 aligner bundled for R
    "Rbowtie2", # Bowtie2 aligner bundled for R
    "QuasR"     # higher-level wrapper: preprocessing + alignment + quantification
)
install_if_missing(alignment_pkgs)

## ---- 3. BAM handling & quantification --------------------------------------

message("\n[Section] BAM handling & quantification ...")
quant_pkgs <- c(
    "Rsamtools",
    "GenomicAlignments",
    "GenomicFeatures",
    "GenomicRanges",
    "SummarizedExperiment",
    "tximport", # import transcript-level quantifications (e.g. external Salmon/kallisto output)
    "tximeta"   # tximport wrapper with automatic metadata linking
)
install_if_missing(quant_pkgs)

## ---- 4. Human annotation & genome resources --------------------------------

message("\n[Section] Human annotation & genome resources ...")
annotation_pkgs <- c(
    "org.Hs.eg.db",                     # gene ID mapping
    "TxDb.Hsapiens.UCSC.hg38.knownGene", # transcript database (hg38)
    "EnsDb.Hsapiens.v86",                # Ensembl-based annotation
    "AnnotationDbi",
    "AnnotationHub",
    "ensembldb",
    "biomaRt"
)
install_if_missing(annotation_pkgs)

## NOTE: BSgenome.Hsapiens.UCSC.hg38 is a large download (~800MB+).
## Comment this block out if you don't need full genome sequence (e.g. you're
## only doing transcript-level quantification with tximport/tximeta).
genome_pkgs <- c(
    "BSgenome.Hsapiens.UCSC.hg38"
)
install_if_missing(genome_pkgs)

## ---- 5. Differential expression (DESeq2 primary) ---------------------------

message("\n[Section] Differential expression ...")
de_pkgs <- c(
    "DESeq2", # primary DE tool
    "apeglm", # effect-size shrinkage for lfcShrink()
    "ashr",   # alternative shrinkage estimator
    "IHW",    # independent hypothesis weighting for multiple testing
    "vsn"     # variance stabilization diagnostics
)
install_if_missing(de_pkgs)

## ---- 6. Downstream analysis / visualization / enrichment -------------------

message("\n[Section] Downstream analysis, visualization & enrichment ...")
downstream_pkgs <- c(
    "pheatmap",       # heatmaps (CRAN)
    "EnhancedVolcano", # volcano plots
    "PCAtools",       # PCA diagnostics
    "clusterProfiler", # GO/KEGG enrichment
    "enrichplot",     # enrichment result visualization
    "DOSE",           # disease ontology enrichment
    "ReactomePA"      # Reactome pathway enrichment
)
install_if_missing(downstream_pkgs)

## ---- 7. General data wrangling / reporting (CRAN) --------------------------

message("\n[Section] Data wrangling & reporting ...")
general_pkgs <- c(
    "tidyverse",
    "rmarkdown",
    "knitr"
)
install_if_missing(general_pkgs)

## ---- 8. Verification: load every package and report PASS/FAIL -------------

message("\n[Section] Verifying installation ...")

all_pkgs <- c(
    qc_trimming_pkgs,
    alignment_pkgs,
    quant_pkgs,
    annotation_pkgs,
    genome_pkgs,
    de_pkgs,
    downstream_pkgs,
    general_pkgs
)

results <- vapply(all_pkgs, function(pkg) {
    ok <- requireNamespace(pkg, quietly = TRUE)
    if (ok) "PASS" else "FAIL"
}, character(1))

summary_df <- data.frame(package = all_pkgs, status = results, row.names = NULL)
print(summary_df, row.names = FALSE)

n_fail <- sum(results == "FAIL")
if (n_fail > 0) {
    warning(sprintf(
        "%d package(s) failed to install/load. Check the messages above for missing system dependencies (compilers, libxml2, libcurl, etc.).",
        n_fail
    ))
} else {
    message("\nAll packages installed and loaded successfully.")
}
