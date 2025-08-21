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
sample_sheet <- read.table(sample_sheet_file, header = TRUE, row.names=NULL, fill = TRUE, quote = "\"")
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
print(resultsNames(dds)[0])
print(resultsNames(dds)[1])
print(resultsNames(dds)[2])
print(resultsNames(dds)[3])
# Output results
for (conditionName in resultsNames(dds)[-1]) {
  write.csv(as.data.frame(results(dds, name = conditionName)), file = file.path(output_dir, paste0(gsub("[^A-Za-z0-9_]", "_", conditionName),"_deseq2_results.csv")))
}

# Loop through all result names except the intercept
for (conditionName in resultsNames(dds)[-1]) {
  
  # Get DE results
  res <- results(dds, name = conditionName)
  
  # Prepare volcano data
  volcano_data <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(
      neg_log10_padj = -log10(padj),
      significance = case_when(
        padj < 0.05 & abs(log2FoldChange) > 1 ~ "Significant",
        TRUE ~ "Not Significant"
      ),
      neg_log10_padj = ifelse(is.infinite(neg_log10_padj), NA, neg_log10_padj)
    )
  
  # Create volcano plot
  volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = neg_log10_padj)) +
    geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    theme_minimal() +
    theme(
      legend.title = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(
      title = paste("Volcano Plot:", conditionName),
      x = "log2(Fold Change)",
      y = "-log10(padj)"
    )
  
  # Clean the condition name for safe filenames
  file_safe_name <- gsub("[^A-Za-z0-9_]", "_", conditionName)
  
  # Save plot
  ggsave(
    filename = file.path(output_dir, paste0(file_safe_name, "_volcano_plot.png")),
    plot = volcano_plot,
    width = 8,
    height = 6,
    bg = "white"
  )
}

# PCA plot

PCA <- plotPCA(vsd)
ggsave(filename = file.path(output_dir, "pca_plot.png"), plot = PCA, width = 8, height = 6)


# Plot expression for target genes listed in sample_sheet

# Check if "target_genes" column exists
if ("target_genes" %in% colnames(sample_sheet)) {
  
  # Gather all genes from all rows
  target_genes <- unlist(strsplit(as.character(sample_sheet$target_genes), ",")) %>%
    trimws() %>%
    unique()

  # Check which of the target genes are in the dataset
  genes_in_data <- target_genes[target_genes %in% rownames(rld)]
  missing_genes <- setdiff(target_genes, genes_in_data)

  if (length(missing_genes) > 0) {
    warning("Some target genes were not found in the data and will be skipped: ", paste(missing_genes, collapse = ", "))
  }

  # Create output folder for gene plots
  gene_plot_dir <- file.path(output_dir, "gene_expression_plots")
  dir.create(gene_plot_dir, showWarnings = FALSE)

  # Plot expression for each gene
  for (gene in genes_in_data) {
    df <- plotCounts(dds, gene = gene, intgroup = "condition", returnData = TRUE)
    
    p <- ggplot(df, aes(x = condition, y = count)) +
      geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
      stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "red") +
      theme_minimal() +
      labs(title = paste("Expression of", gene),
           y = "Normalized count",
           x = "Condition")
    
    ggsave(filename = file.path(gene_plot_dir, paste0(gene, "_expression.png")),
           plot = p, width = 6, height = 4, bg="white")
  }
}

# Plot expression using raw counts for target genes

if ("target_genes" %in% colnames(sample_sheet)) {
  
  # Extract all unique target genes from the sample sheet
  target_genes <- unlist(strsplit(as.character(sample_sheet$target_genes), ",")) %>%
    trimws() %>%
    unique()

  # Filter genes that are present in the counts matrix
  genes_in_data <- target_genes[target_genes %in% rownames(counts)]
  missing_genes <- setdiff(target_genes, genes_in_data)

  if (length(missing_genes) > 0) {
    warning("Some target genes were not found in the count matrix: ", paste(missing_genes, collapse = ", "))
  }

  # Create output folder for gene expression plots
  gene_plot_dir <- file.path(output_dir, "gene_expression_plots_raw")
  dir.create(gene_plot_dir, showWarnings = FALSE)

  # Loop over genes and create a plot for each
  for (gene in genes_in_data) {
  # Get the raw counts for the gene and reshape
    gene_counts <- data.frame(
    sample = colnames(counts),
    count = as.numeric(counts[gene, ])
    )

  # Join with sample_sheet to get condition
    gene_data <- left_join(gene_counts, sample_sheet, by = "sample")

    p <- ggplot(gene_data, aes(x = condition, y = count)) +
     geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
     stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "red") +
     theme_minimal() +
     scale_y_log10() +
     labs(title = paste("Raw Counts for", gene),
          y = "Raw count",
          x = "Condition")

  # Save plot
    ggsave(filename = file.path(gene_plot_dir, paste0(gene, "_raw_expression.png")),
         plot = p, width = 6, height = 4, bg="white")
  }

}
