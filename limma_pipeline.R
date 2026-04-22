#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(limma); library(edgeR); library(sva); library(clusterProfiler); library(org.Hs.eg.db); library(tidyverse); library(EnhancedVolcano)
})
cat('=== limma-voom Differential Expression Pipeline ===\n')
counts <- read.csv('data/counts.csv', row.names = 1, check.names = FALSE)
design_df <- read.csv('data/design.csv', row.names = 1)
dge <- DGEList(counts = counts)
keep <- filterByExpr(dge, group = design_df$condition)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = 'TMM')
if ('batch' %in% colnames(design_df)) {
  cat('Applying ComBat-seq batch correction...\n')
  dge$counts <- ComBat_seq(counts = as.matrix(dge$counts), batch = design_df$batch, group = design_df$condition)
}
design_matrix <- model.matrix(~ condition, data = design_df)
v <- voom(dge, design_matrix, plot = FALSE)
fit <- lmFit(v, design_matrix)
fit <- eBayes(fit)
results <- topTable(fit, coef = 2, number = Inf, sort.by = 'P')
sig <- results %>% filter(adj.P.Val < 0.05)
write.csv(sig, 'results/limma_DEGs.csv')
cat(sprintf('Significant DEGs: %d\n', nrow(sig)))
