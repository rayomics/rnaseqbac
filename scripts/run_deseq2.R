#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
count_file <- args[1]
sample_sheet_file <- args[2]
output_dir <- args[3]

suppressMessages({
  library(DESeq2)
  library(tidyverse)
})

# Read count matrix
counts <- read.delim(count_file, comment.char = "#", row.names = 1)
# Keep only sample columns
counts <- counts[, 6:ncol(counts)]

# Extract sample names from BAM paths 
colnames(counts) <- sub(".*\\.(.*)\\.bam$", "\\1", colnames(counts))

# Read sample sheet
sample_sheet <- read.table(sample_sheet_file, header = TRUE, row.names=NULL)
sample_sheet$sample <- sub("_R[12]_.*", "", sample_sheet$file_R1)

# Check if all samples in sample_sheet are in the counts matrix
if (!all(sample_sheet$sample %in% colnames(counts))) {
  missing <- setdiff(sample_sheet$sample, colnames(counts))
  stop("The following samples are missing from the count matrix: ", paste(missing, collapse = ", "))
}

# Subset and reorder counts to match sample_sheet
counts <- counts[, sample_sheet$sample]
head(counts)
# Create coldata
coldata <- data.frame(row.names = sample_sheet$sample,
                      condition = sample_sheet$condition)
print(sample_sheet)
# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

# Output results
write.csv(as.data.frame(res), file = file.path(output_dir, "deseq2_results.csv"))
