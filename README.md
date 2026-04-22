# limma-voom Analysis Pipeline

> You've got 3 samples per group. DESeq2 starts to wobble with sample sizes this small. limma doesn't. Here's why, and how to use it.

## The Problem This Solves

Not every experiment has the luxury of 10 biological replicates. Clinical samples are expensive. Patient cohorts are limited. limma, with its empirical Bayes moderation, borrows information across genes to stabilise variance estimates, making it the better choice when your sample size is tight.

This pipeline also handles batch effects. You ran half your samples in January and half in March? ComBat-seq handles this before it ruins your DEGs.

## What Happens When You Run This

1. **TMM normalisation** corrects for compositional differences between libraries.
2. **ComBat-seq batch correction** removes technical variation while preserving biological signal. This is the step most people skip and then wonder why their PCA clusters by sequencing batch instead of condition.
3. **voom transformation** converts counts to log-CPM with precision weights, letting limma work beautifully with RNA-seq data.
4. **Empirical Bayes moderation** shrinks gene-wise variances toward a common value, which is exactly what you need when n=3.
5. **Gene set enrichment** with clusterProfiler gives you the biological story, not just a list of genes.

## Quick Start
```bash
Rscript limma_pipeline.R --input data/expression_matrix.csv --design data/design.csv
```

## A Real Decision Framework

| Scenario | Use This |
|----------|----------|
| n >= 5 per group, RNA-seq | DESeq2 or limma-voom (both fine) |
| n < 5 per group | **limma-voom** (this pipeline) |
| Batch effects present | **This pipeline** (ComBat-seq built in) |
| Complex design (interaction terms) | limma handles this more naturally than DESeq2 |

## The Uncomfortable Truth About Batch Effects

If you don't correct for batch effects, your "top differentially expressed gene" might just be the gene most affected by which day you ran your experiment. Run the PCA first, colour by batch, and if the samples cluster by batch before condition, you need ComBat before anything else.
