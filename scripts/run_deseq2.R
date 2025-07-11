#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
count_file <- args[1]
output_dir <- args[2]

suppressMessages({
  library(DESeq2)
  library(tidyverse)
})

counts <- read.delim(count_file, comment.char = "#", row.names = 1)
counts <- counts[,6:ncol(counts)]  # remove annotation columns
coldata <- data.frame(row.names=colnames(counts),
                      condition=rep(c("control", "treated"), each=3))  # adjust as needed

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

write.csv(as.data.frame(res), file = file.path(output_dir, "deseq2_results.csv"))
