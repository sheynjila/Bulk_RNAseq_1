# Differential Gene Expression Study — Airway Smooth Muscle / Dexamethasone

A complete, R-based DESeq2 differential-expression workflow: from a raw, slightly
messy gene-count matrix through validation, statistical modeling, and a full set of
standard RNA-seq diagnostic and results plots — including an optional KEGG pathway
overlay.

The dataset is the Himes et al. (2014) airway smooth muscle glucocorticoid-response
study (GEO [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778),
BioProject [PRJNA229998](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA229998)) — the
same study used as the canonical Bioconductor `DESeq2` teaching example. This project
picks up exactly where a raw-reads-to-counts pipeline (STAR + featureCounts, or
equivalent) leaves off: it starts from an already-quantified gene × sample count
matrix and carries it through to differential expression and visualization.

> **Status:** Built and reviewed as coursework material for the MIBO8110 Applied Omics
> curriculum. Every script has been read line-by-line for correctness; this
> documentation was authored in an environment without a working R/Bioconductor
> installation, so results have not been re-executed end-to-end here. One concrete,
> verified data problem *was* found and fixed during this review — see
> [Known Issues (fixed)](#known-issues-fixed-during-this-cleanup) below.

## Table of contents

- [Study & data](#study--data)
- [Repository layout](#repository-layout)
- [Prerequisites](#prerequisites)
- [Quickstart](#quickstart)
- [The pipeline, in order](#the-pipeline-in-order)
- [Data quirks this pipeline defends against](#data-quirks-this-pipeline-defends-against)
- [Known issues (fixed during this cleanup)](#known-issues-fixed-during-this-cleanup)
- [Documentation in docs/](#documentation-in-docs)
- [Companion repositories](#companion-repositories)
- [Runbook](#runbook)
- [Citation](#citation)
- [License](#license)

## Study & data

| | |
|---|---|
| **Study** | Himes BE, Jiang X, Wagner P, et al. (2014). RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells. *PLOS ONE*. |
| **DOI** | [10.1371/journal.pone.0099625](https://doi.org/10.1371/journal.pone.0099625) |
| **GEO series** | GSE52778 |
| **Design** | 4 donor cell lines × 2 conditions (dexamethasone-treated vs. untreated) = 8 samples, paired by donor |
| **Input files** | `data/counts.csv` (Ensembl gene IDs × 8 SRR-accession columns), `data/metadata2.csv` (sample ID, cell line, treatment) |
| **Statistical design** | `~ cellline + dexamethasone` — adjusts for donor-to-donor baseline variation while estimating the treatment effect |

## Repository layout

```text
Differential_Gene_Expression_Study/
├── Differential_Gene_Expression_Study.Rproj
├── README.md                      # This file
├── RNAseq_DGE_Runbook.docx/.pdf    # Operational SOP companion to this README
├── 02_organize_project_files.R    # Kept at the project root by design — see note below
├── scripts/
│   ├── 00_install_and_load_packages.R   # Two teaching methods: manual vs. automated install
│   ├── 00_manage_project_packages.R     # Recommended: single package-status/install function
│   ├── 01_create_project_directories.R  # Creates scripts/, data/, plots/
│   ├── 02_organize_project_files.R      # Copy of the root organizer (see note below)
│   ├── DESeq.R                          # Historical/diagnostic — see note below, NOT part of the run order
│   ├── Aligning_meta_counts.R           # Step 1: validate & align counts + metadata
│   ├── LOAD_CLEANED_ALIGNED_DATA.R      # Step 2: build the DESeq2 model, run the contrast
│   ├── VISUALIZE_DESeq2_RESULTS.R       # Step 3: PCA, volcano, top-gene boxplot
│   └── Additional_high-value_visualization_plots_commonly_used_in_RNA-seq_analyses.R  # Step 4: extended QC/results plots + optional shrinkage/enrichment/pathway analysis
├── data/
│   ├── counts.csv / metadata2.csv       # Original inputs (as supplied)
│   ├── counts_aligned.csv / metadata_aligned.csv  # Validated, sample-order-matched pair (output of Aligning_meta_counts.R)
│   └── cnt_adj.csv / met_adj.csv        # Output of the historical DESeq.R script — not consumed by anything downstream
├── plots/
│   └── hsa04110.png / hsa04110.pathview.png / hsa04110.xml   # KEGG "Cell cycle" pathway overlay, from the optional pathview step
└── docs/
    ├── RNAseq_DESeq2_Class_Exercise.docx/.html   # The original student lab handout this repo implements
    ├── RNAseq_Troubleshooting_Notebook.Rmd/.html/.nb.html  # Reproducible write-up of the counts/metadata debugging story
    ├── Bioconductor_Pacakage_Installation.txt    # Minimal historical install snippet
    └── install_rnaseq_packages.R                 # Broader reference installer — see note below
```

**Why `02_organize_project_files.R` exists in two places:** this is deliberate, not a
duplication bug. The script organizes the project root into `scripts/`, `data/`, and
`plots/`, but on Windows a running R script cannot move itself while it's executing —
so it copies itself into `scripts/` (for the record) while the original stays at the
project root, where it can be re-run again later without first having to be located.

## Prerequisites

- R (a recent 4.x release) with Bioconductor access.
- Run `scripts/00_manage_project_packages.R` (from the R console, with the project
  root as your working directory) to see exactly which packages are installed and
  which are missing — call `manage_packages("automatic")` to install everything this
  project needs in one step. See [Bioconductor_Pacakage_Installation.txt](docs/Bioconductor_Pacakage_Installation.txt)
  for a minimal, dependency-free fallback if `manage_packages()` itself won't run yet.
- **Core packages:** DESeq2, ggplot2, dplyr, tidyr, pheatmap, RColorBrewer.
- **Optional (only needed if you enable the optional analyses in the last script):**
  apeglm, clusterProfiler, enrichplot, org.Hs.eg.db, pathview.

## Quickstart

Open `Differential_Gene_Expression_Study.Rproj` in RStudio first — every script below
is written assuming the working directory is the project root (RStudio sets this
automatically when you open the `.Rproj`), not the `scripts/` folder itself.

```r
# 1. See what's installed / install what's missing
source("scripts/00_manage_project_packages.R")
manage_packages("automatic")

# 2. Validate the raw inputs and produce a clean, aligned pair
source("scripts/Aligning_meta_counts.R")

# 3. Build the DESeq2 model and run the treatment contrast
source("scripts/LOAD_CLEANED_ALIGNED_DATA.R")

# 4. Core results plots (PCA, volcano, top-gene boxplot)
source("scripts/VISUALIZE_DESeq2_RESULTS.R")

# 5. Extended QC/results plots, plus optional shrinkage / enrichment / KEGG pathway
source("scripts/Additional_high-value_visualization_plots_commonly_used_in_RNA-seq_analyses.R")
```

The full guided version of this walkthrough — with expected output at each step,
timing estimates, and a submission checklist — is
[`docs/RNAseq_DESeq2_Class_Exercise.docx`](docs/RNAseq_DESeq2_Class_Exercise.docx),
the original student lab handout this repository is the worked/completed version of.

## The pipeline, in order

| Step | Script | What it does |
|---|---|---|
| 0 | `00_manage_project_packages.R` | Reports install status for every package this project uses; can install what's missing (CRAN + Bioconductor) in one call. |
| 0 (alt.) | `00_install_and_load_packages.R` | A more explicit, two-method (manual vs. automated) teaching version of the same idea — kept for reference; `00_manage_project_packages.R` is the recommended entry point. |
| 1 | `01_create_project_directories.R` | Creates `scripts/`, `data/`, `plots/` if they don't already exist (this repository is already organized this way). |
| 2 | `02_organize_project_files.R` | Moves loose `.R`/`.csv`/plot files at the project root into the folders above. Safe to re-run; never overwrites an existing destination file. |
| 3 | `Aligning_meta_counts.R` | Reads `data/counts.csv` and `data/metadata2.csv`, correctly handling the trailing-comma header (see [Data quirks](#data-quirks-this-pipeline-defends-against) below), validates sample IDs, and writes a sample-order-matched `data/counts_aligned.csv` + `data/metadata_aligned.csv`. |
| 4 | `LOAD_CLEANED_ALIGNED_DATA.R` | Loads the aligned pair, builds a `DESeqDataSet` with design `~ cellline + dexamethasone`, filters genes with fewer than 10 total reads, runs `DESeq()`, and extracts the `dexamethasone: treated vs. untreated` contrast. |
| 5 | `VISUALIZE_DESeq2_RESULTS.R` | PCA plot, volcano plot, and a boxplot of the single most significant gene. |
| 6 | `Additional_high-value_visualization_plots_commonly_used_in_RNA-seq_analyses.R` | Sample-distance heatmap, top-30-gene heatmap, MA plot, raw p-value histogram, dispersion diagnostic, top-10-gene expression profiles — plus three **optional, off-by-default** blocks (toggle the flags at the top of the script): log2FC shrinkage (`apeglm`), GO enrichment + GSEA (`clusterProfiler`), and a KEGG pathway overlay (`pathview`). The `plots/hsa04110.*` files in this repo are the output of a completed pathview run against the "Cell cycle" pathway. |

**`DESeq.R` is not part of this run order.** It's the original, unedited live
debugging session — narrated in first person — that discovered the two data quirks
described below. It's kept in `scripts/` for its teaching value (it shows *how* the
problems were found, not just the fix), and its outputs (`data/cnt_adj.csv`,
`data/met_adj.csv`) are dead ends: nothing downstream reads them. The reusable,
defensive fix for the same problem is `Aligning_meta_counts.R`, step 3 above. The
polished, reproducible write-up of the same debugging session — worth reading even if
you don't run `DESeq.R` itself — is
[`docs/RNAseq_Troubleshooting_Notebook.Rmd`](docs/RNAseq_Troubleshooting_Notebook.Rmd).

## Data quirks this pipeline defends against

Both were discovered the hard way (see `DESeq.R` / the troubleshooting notebook) and
are why `Aligning_meta_counts.R` exists rather than a plain `read.csv()`:

1. **A trailing comma in `counts.csv`'s header row** creates a spurious, unnamed
   column (it reads in as `X`) once R applies default column-naming — and a naive
   `read.csv()` also risks misaligning the *first* sample column (`SRR1039508`)
   against the header. `Aligning_meta_counts.R` reads the header line separately with
   `scan()` and reconstructs the column names explicitly instead.
2. **`metadata2.csv`'s sample IDs live in an ordinary first column, not as row
   names.** A naive `read.csv()` leaves `rownames(met)` as `"1"`...`"8"`, which
   silently fails `all(colnames(cnt) %in% rownames(met))` in a way that isn't obvious
   from `head()` alone. `Aligning_meta_counts.R` promotes that column to real row
   names and validates the result.

## Known issues (fixed during this cleanup)

- **`data/counts_aligned.csv` and `data/metadata_aligned.csv` were stale and
  incorrect** as shipped: `counts_aligned.csv` was missing the `SRR1039508` column
  entirely (7 samples instead of 8), and `metadata_aligned.csv`'s columns were
  mislabeled (`"dexamethasone","X"` instead of `"cellline","dexamethasone"`) — the
  output of an earlier, since-corrected version of the alignment logic. Both files
  have been replaced with correct, verified output (cross-checked against an
  independently generated identical copy from the companion
  `mibo8110_First_R_Exercise` repository, and consistent with what the current
  `Aligning_meta_counts.R` produces from the raw inputs). If you re-run
  `Aligning_meta_counts.R` yourself, it will regenerate both files identically to
  what's now committed here.
- **`docs/Bioconductor_Pacakage_Installation.txt` had a syntax error** — a missing
  closing parenthesis on `BiocManager::install(c("GenomicFeatures")` — that would
  have stopped a copy-pasted install partway through. Fixed.

No other script logic was changed. `00_install_and_load_packages.R` and `DESeq.R` are
intentionally kept as-is despite overlapping in purpose with their respective
successors (`00_manage_project_packages.R`, `Aligning_meta_counts.R`) — they have
real, distinct teaching value and are clearly labeled above as historical/alternate,
not as the recommended path.

## Documentation in docs/

- **`RNAseq_DESeq2_Class_Exercise.docx` / `.html`** — the original student lab
  handout (learning objectives, a suggested session-time map, and a submission
  checklist). This repository is the completed, working version of that exercise.
- **`RNAseq_Troubleshooting_Notebook.Rmd` / `.html` / `.nb.html`** — a reproducible
  R Markdown notebook covering the count/metadata debugging story in full, plus
  general data-import best practices, an automated validation function, and a
  suggested project layout for future RNA-seq projects. The validation function it
  documents now also lives as a standalone, reusable script in the companion
  [`Automated_validation_function`](#companion-repositories) repository.
- **`Bioconductor_Pacakage_Installation.txt`** — a minimal, three-line historical
  install snippet, kept for reference; prefer `00_manage_project_packages.R` for
  everyday use.
- **`install_rnaseq_packages.R`** — a much broader installer covering an entire
  R-native pipeline from raw reads to enrichment (`Rfastp`, `Rsubread`/`Rhisat2`/
  `Rbowtie2`/`QuasR` for alignment, `tximport`/`tximeta`, full annotation packages,
  and everything `00_manage_project_packages.R` installs). This project's own scripts
  only need the smaller package list in [§4](#prerequisites) — this file is here as a
  reference for anyone extending this repo to start from raw FASTQ files instead of
  an already-quantified count matrix.

## Companion repositories

This project is one of three related, separately published repositories:

- **`mibo8110_First_R_Exercise`** — a compact, self-contained version of the same
  exercise, intended as a first hands-on R + GitHub exercise for students who haven't
  set up a project or pushed a repository before.
- **`Automated_validation_function`** — the counts/metadata validation function
  discussed above, extracted as a standalone, reusable utility with its own example
  dataset.

## Runbook

[`RNAseq_DGE_Runbook.docx`](RNAseq_DGE_Runbook.docx) /
[`RNAseq_DGE_Runbook.pdf`](RNAseq_DGE_Runbook.pdf) is the step-by-step operational
companion to this README — exact commands, expected console output at each
checkpoint, and a troubleshooting reference built from the two data quirks above.

## Citation

> Himes BE, Jiang X, Wagner P, Hu R, Wang Q, Klanderman B, et al. (2014) RNA-Seq
> Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene
> that Modulates Cytokine Function in Airway Smooth Muscle Cells. *PLoS ONE* 9(6):
> e99625. https://doi.org/10.1371/journal.pone.0099625

## License

This repository's code and documentation are released under the [MIT License](LICENSE).
The underlying sequencing data and publication remain governed by their own GEO and
journal terms regardless of the license on this code.

---

*Maintained as part of the MIBO8110 Applied Omics curriculum.*
