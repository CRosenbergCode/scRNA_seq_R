#https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/slingshot/slingshot.html

BiocParallel::register(BiocParallel::SerialParam())
#BiocManager::install("slingshot")


library(slingshot)
#library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)

seu = readRDS("rdsFiles/seurat_clustered_1_23.RDS")


dimred <- seu@reductions$umap@cell.embeddings
clustering <- seu$RNA_snn_res.0.8
counts <- as.matrix(seu@assays$RNA@counts[seu@assays$RNA@var.features, ])

set.seed(1)
lineages <- getLineages(data = dimred, clusterLabels = clustering,start.clus = '0')

lineages

pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))


par(mfrow = c(1, 2))
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
for (i in levels(clustering)) {
  text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
}
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
lines(SlingshotDataSet(lineages), lwd = 3, col = "black")

curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves

plot(dimred, col = pal[clustering], asp = 1, pch = 16)
lines(SlingshotDataSet(curves), lwd = 3, col = "black")



BiocManager::install("tradeSeq",version='3.19')
#BiocManager::install(version = '3.19')
library(tradeSeq)

# Removing some genes to speed up the computations for this tutorial
filt_counts <- counts[rowSums(counts > 5) > ncol(counts)/100, ]
dim(filt_counts)


sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves)

plotGeneCount(curves, filt_counts, clusters = clustering, models = sce)


#counts <- as.matrix(seu@assays$RNA@counts[seu@assays$RNA@var.features, ])

counts <- as.matrix(seu@assays$RNA@counts)


#counts <- as.matrix(seu[["RNA"]]$counts)

filt_counts <- counts[rowSums(counts > 5) > ncol(counts)/100, ]

saveRDS(filt_counts,file="GAMfitFILTCounts.RDS")


#sce <- fitGAM(counts = as.matrix(counts), sds = curves,nknots=15)
#sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves,nknots=15)

plotGeneCount(curves, filt_counts, clusters = clustering, models = sce)

#saveRDS(sce,file="GAMfitSCE.RDS")
sce=readRDS(file="GAMfitSCE.RDS")
pseudotime_association <- associationTest(sce)
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)

library(dplyr)
plot_differential_expression <- function(feature_id) {
  feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  cowplot::plot_grid(plotGeneCount(curves, filt_counts, gene = feature_id[1], clusters = clustering, models = sce) + ggplot2::theme(legend.position = "none"), 
                     plotSmoothers(sce, as.matrix(counts), gene = feature_id[1]))
}

feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
plot_differential_expression(feature_id)
