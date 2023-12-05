library(Seurat)
library(dplyr)
library(singleCellTK)
library(dittoSeq)
#library(SeuratDisk)
#library(SeuratData)
#packageDescription("singleCellTK")$Version
#packageDescription("SummarizedExperiment")$Version
#R.version
#singleCellTK()
cxtovary = Read10X(data.dir = "~/datasets/sample1")
cxt_seurat = CreateSeuratObject(counts = cxtovary)
ctx_seurat = NormalizeData(object = cxt_seurat)
cxt_seurat = FindVariableFeatures(object = cxt_seurat)
cxt_seurat  = ScaleData(cxt_seurat, features = rownames(var_feats)) #Remove features = to increase speed but reduce number of scaled genes
cxt_seurat  = RunPCA(cxt_seurat )
cxt_seurat  = FindNeighbors(cxt_seurat )
cxt_seurat  = FindClusters(object = cxt_seurat )
cxt_seurat  = RunTSNE(cxt_seurat )

cxt_single = as.SingleCellExperiment(tsnes)
cxt_non_proc = as.SingleCellExperiment(cxt_seurat)
sce = runCellQC(cxt_non_proc, sample = NULL,
                algorithms = c("QCMetrics", "scDblFinder", "decontX"),
                seed = 12345)#,mitoID=mtGenes,mitoIDType = "symbol")

umap_sce = runQuickUMAP(sce, reducedDimName = "QC_UMAP",seed = 2023, sample = NULL)

#sce_cols = subsetSCECols(sce, colData = c("total > 75", 
#                                      "detected > 20"))
sce_cols = subsetSCECols(sce, colData = c("total > 600", 
                                          "detected > 300"))

seur_sce = runSeuratNormalizeData(inSCE = sce, useAssay = "decontXcounts", normAssayName = "seuratNormData", normalizationMethod = "LogNormalize", scaleFactor = 10000)

seur_sce = runSeuratFindHVG(inSCE = seur_sce, useAssay = "decontXcounts", method = "vst", hvgNumber = 2000, createFeatureSubset = "hvf")

# Print names of top 10 variable features
print(getTopHVG(inSCE = seur_sce, method = "vst", hvgNumber = 10))

seur_sce = runSeuratPCA(inSCE = seur_sce, useAssay = "seuratNormData", reducedDimName = "pca", nPCs = 50, seed = 42, scale = TRUE, useFeatureSubset = "hvf")

#seur_sce <- runSeuratJackStraw(inSCE = seur_sce, useAssay = "seuratNormData", dims = 50)

# Plot JackStraw
#plotSeuratJackStraw(inSCE = seur_sce, dims = 50)

seur_sce = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 10) 
# sce_2 <- findMarkerDiffExp(seur_sce, useAssay = "logcounts", method = "wilcox",
#                          cluster = "cluster",
#                          log2fcThreshold = 0, fdrThreshold = 0.05,
#                          minClustExprPerc = 0, maxCtrlExprPerc = 1,
#                          minMeanExpr = 0)
# markerHm <- plotFindMarkerHeatmap(seur_sce, topN = 5, log2fcThreshold = 0, 
#                               fdrThreshold = 0.05, minClustExprPerc = 0.5, 
#                               maxCtrlExprPerc = 0.4, minMeanExpr = 0, 
#                               rowLabel = TRUE)'''

#sub_clustering()

seur_sce = runSeuratFindMarkers(inSCE = seur_sce, allGroup = "Seurat_louvain_Resolution0.8")

# Fetch marker genes table
markerGenes = metadata(seur_sce)[["seuratMarkers"]] 

# Order by log fold change and p value
markerGenes = markerGenes[order(-markerGenes$avg_log2FC, markerGenes$p_val),]

markerGenes1 = markerGenes

head(markerGenes)
plotSeuratGenes(inSCE = seur_sce, useAssay = "seuratNormData", plotType = "ridge", features = metadata(seur_sce)[["seuratMarkers"]]$gene.id[9:12], groupVariable = "Seurat_louvain_Resolution0.8", ncol = 2, combine = TRUE)
FeaturePlot(seur_sce, features = features)

#cxtovary = Read10X(data.dir = "~/datasets/sample3")
#cxt_seurat = CreateSeuratObject(counts = cxtovary)

#colnames(cxt_seurat)

temp_seur = as.Seurat(seur_sce, counts = "counts", data = "logcounts")

temp_seur <- FindNeighbors(temp_seur, dims = 1:10)
temp_seur <- FindClusters(temp_seur)

temp_seur = clusts

saveRDS(temp_seur, file = "sample3clustered.rds")

temp_seur.markers <- FindAllMarkers(temp_seur, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


temp_seur.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(temp_seur, features = top10$gene) + NoLegend()

head(Idents(temp_seur), 5)


DoHeatmap(object = seur_sce)

dittoHeatmap(seur_sce, genes,
             annot.by = "clustering")

genes <- getGenes(seur_sce)[1:100]

plotSCEHeatmap(seur_sce)

sce <- runScranSNN(sce, useReducedDim = "PCA", clusterName = "cluster", algorithm = "louvain", k = 4)

sce <- findMarkerDiffExp(sce, useAssay = "logcounts", method = "wilcox",
                         cluster = "cluster",
                         log2fcThreshold = 0, fdrThreshold = 0.05,
                         minClustExprPerc = 0, maxCtrlExprPerc = 1,
                         minMeanExpr = 0)

markerHm <- plotMarkerDiffExp(seur_sce, topN = 5, log2fcThreshold = 0, 
                              fdrThreshold = 0.05, minClustExprPerc = 0.5, 
                              maxCtrlExprPerc = 0.4, minMeanExpr = 0, 
                              rowLabel = TRUE)

markerHm