#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
count_file <- args[1]
sample_sheet_file <- args[2]
output_dir <- args[3]

# Ensure GenomeInfoDbData is installed
if (!requireNamespace("GenomeInfoDbData", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install("GenomeInfoDbData", ask = FALSE, update = FALSE)
}

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
#sample_sheet$sample <- sub("_R[12]_.*", "", sample_sheet$file_R1)
sample_sheet$sample <- sub("_R1.*", "", sample_sheet$file_R1)

# Check if all samples in sample_sheet are in the counts matrix
if (!all(sample_sheet$sample %in% colnames(counts))) {
  missing <- setdiff(sample_sheet$sample, colnames(counts))
  stop("The following samples are missing from the count matrix: ", paste(missing, collapse = ", "))
}

# Subset and reorder counts to match sample_sheet
counts <- counts[, sample_sheet$sample]

# Create coldata
coldata <- data.frame(row.names = sample_sheet$sample,
                      condition = sample_sheet$condition)

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

# Output results
write.csv(as.data.frame(res), file = file.path(output_dir, "deseq2_results.csv"))

# Create volcano plot
volcano_data <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  mutate(
    neg_log10_padj = -log10(padj),
    significance = case_when(
      padj < 0.05 & abs(log2FoldChange) > 1 ~ "Significant",
      TRUE ~ "Not Significant"
    )
  )

# Replace Inf/-Inf with NA
volcano_data <- volcano_data %>%
  mutate(
    neg_log10_padj = ifelse(is.infinite(neg_log10_padj), NA, neg_log10_padj)
  )

# Plot
volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "log2(Fold Change)", y = "-log10(padj)") +
  theme(legend.title = element_blank())

# Save plot
ggsave(filename = file.path(output_dir, "volcano_plot.png"), plot = volcano_plot, width = 8, height = 6)

# PCA plot

PCA <- plotPCA(vsd)
ggsave(filename = file.path(output_dir, "pca_plot.png"), plot = PCA, width = 8, height = 6)
