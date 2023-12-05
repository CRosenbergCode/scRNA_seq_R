library(Seurat)
library(dplyr)
library(singleCellTK)
sce = readRDS("sample1postqc.rds")

ncol(sce)


#sce_cols = subsetSCECols(sce, colData = c("total > 75", 
#                                          "detected > 20"))
sce_cols = subsetSCECols(sce, colData = c("total > 600", 
                                          "detected > 300"))
#Change to umap_sce when working
ncol(sce_cols)

sce <- runNormalization(sce, useAssay = "counts", outAssayName = "logcounts",
                        normalizationMethod = "logNormCounts")
sce <- runNormalization(sce, useAssay = "logcounts",
                        outAssayName = "logcounts_scaled", scale = TRUE)
sce <- runSeuratFindHVG(inSCE = seur_sce, useAssay = "decontXcounts", method = "vst", hvgNumber = 2000, createFeatureSubset = "hvf")
#sce <- getTopHVG(inSCE = seur_sce, method = "vst", hvgNumber = 100)
plotTopHVG(sce, method = "vst")
sce <- runDimReduce(inSCE = sce, method = "scaterPCA", 
                    useAssay = "logcounts_scaled", useAltExp = "hvg",
                    reducedDimName = "PCA", seed = 12345)
plotDimRed(sce, useReduction = "PCA")
sce <- runScranSNN(sce, useReducedDim = "PCA", clusterName = "cluster", algorithm = "louvain", k = 4)

sce <- runDimReduce(inSCE = sce, method = "scaterUMAP", useReducedDim = "PCA", reducedDimName = "UMAP", seed = 12345)
plotSCEDimReduceColData(sce, colorBy = "cluster", reducedDimName = "UMAP")
sce <- findMarkerDiffExp(sce, useAssay = "logcounts", method = "wilcox",
                         cluster = "cluster",
                         log2fcThreshold = 0, fdrThreshold = 0.05,
                         minClustExprPerc = 0, maxCtrlExprPerc = 1,
                         minMeanExpr = 0)
markerHm <- plotMarkerDiffExp(sce, topN = 5, log2fcThreshold = 0, 
                              fdrThreshold = 0.05, minClustExprPerc = 0.5, 
                              maxCtrlExprPerc = 0.4, minMeanExpr = 0, 
                              rowLabel = TRUE)