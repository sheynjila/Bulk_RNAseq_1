# Installation of Bioconductor (https://www.bioconductor.org/install/) Installers and Packages 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install() 
# Optional: BiocManager::install(version = "3.23")

BiocManager::install(c("GenomicFeatures"))

#Installation of DESeq2                     
BiocManager::install("DESeq2")
getwd()

cnt <- read.csv("counts.csv")
str(cnt)
met <- read.csv("metadata2.csv")
str(met)

# Making sure that row names in the colData matches to column names in the counts_data
all(colnames(cnt) %in% rownames(met))

# Try these diagnostics:
# Samples in counts but not metadata
setdiff(colnames(cnt), rownames(met))  # Output:  "SRR1039508" "SRR1039509" "SRR1039512" "SRR1039513" "SRR1039516" "SRR1039517" "SRR1039520" "SRR1039521" "X" 
#I got  a FALSE feedback That means at least one sample name in your count matrix (cnt) is not present in the metadata row names (met).

# Samples in metadata but not counts
setdiff(rownames(met), colnames(cnt)) # Output: [1] "1" "2" "3" "4" "5" "6" "7" "8"
head(met)
head(cnt)
dim(cnt)
colnames(cnt)

# My  count matrix has 9 columns with an unwanted x column whiles my metadata has only the 8 real samples, which is why
# all(colnames(cnt) %in% rownames(met)) returns FALSE

# FIX: Remove the spurious X column:
cnt <- cnt[, colnames(cnt) != "X"]
# Check:
dim(cnt)

#REDO: Making sure that row names in the colData matches to column names in the counts_data
all(colnames(cnt) %in% rownames(met))
# Since removing X didn't fix it, let's inspect the row names in met directly.
rownames(met)
head(rownames(met))
# The sample IDs are currently stored in the first column of met, not in the row names.
# Check:
colnames(met)
# Then assign the first column as row names:
rownames(met) <- met[,1]
met <- met[,-1]

# run head(met)
head(rownames(met))

# Now test again:
all(colnames(cnt) %in% rownames(met))
# TRUE
str(met)
str(cnt)


# Want to serve adjusted dataframe in csv

# Install the readr package 
install.packages("readr") 
library(readr)
write_csv(met, "met_adj.csv")
write_csv(cnt, "cnt_adj.csv")
