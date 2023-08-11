library(Seurat)
library(dplyr)
library(singleCellTK)
#library(SeuratDisk)
#library(SeuratData)

#singleCellTK()
cxtovary = Read10X(data.dir = "~/datasets/sample1")
cxt_seurat = CreateSeuratObject(counts = cxtovary)
norm_seurat = NormalizeData(object = cxt_seurat)
var_feats = FindVariableFeatures(object = norm_seurat)
scales = ScaleData(var_feats)
pca = RunPCA(scales)
neighs = FindNeighbors(pca)
clusts = FindClusters(object = neighs)
tsnes = RunTSNE(clusts)
#SaveH5Seurat(tsnes, overwrite = TRUE)
DimPlot(tsnes, reduction = "tsne")

test_single = as.SingleCellExperiment(tsnes)
sce = runCellQC(test_single, sample = NULL,
                 algorithms = c("QCMetrics", "scDblFinder", "decontX"),
                 seed = 12345)

#umap_sce = runQuickUMAP(sce, reducedDimName = "QC_UMAP", seed = 2023)
#Not working yet
#Error in .manageCellVar(inSCE, sample) : 
#Specified variable 'sample' not found in colData(inSCE)

plotRunPerCellQCResults(sce)
#plotScDblFinderResults(umap_sce, reducedDimName = "QC_UMAP")
#Needs umap above to function

#reportCellQC(sce)

hist(sce$nCount_RNA)
hist(sce$nFeature_RNA)

#saveRDS(sce, file = "sample1postqc.rds")
sce = readRDS("sample1postqc.rds")

ncol(sce)


sce_cols = subsetSCECols(sce, colData = c("total > 75", 
                                      "detected > 20"))
#sce_cols = subsetSCECols(sce, colData = c("total > 600", 
#                                          "detected > 300"))
#Change to umap_sce when working
ncol(sce_cols)

seur_sce = runSeuratNormalizeData(inSCE = sce, useAssay = "decontXcounts", normAssayName = "seuratNormData", normalizationMethod = "LogNormalize", scaleFactor = 10000)

seur_sce = runSeuratFindHVG(inSCE = seur_sce, useAssay = "decontXcounts", method = "vst", hvgNumber = 2000, createFeatureSubset = "hvf")

# Print names of top 10 variable features
print(getTopHVG(inSCE = seur_sce, method = "vst", hvgNumber = 10))

plotSeuratHVG(seur_sce, labelPoints = 10)

seur_sce = runSeuratPCA(inSCE = seur_sce, useAssay = "seuratNormData", reducedDimName = "pca", nPCs = 50, seed = 42, scale = TRUE, useFeatureSubset = "hvf")

plotSeuratElbow(inSCE = seur_sce)

seur_sce <- runSeuratJackStraw(inSCE = seur_sce, useAssay = "seuratNormData", dims = 50)

# Plot JackStraw
plotSeuratJackStraw(inSCE = seur_sce, dims = 50)

seur_sce = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 10) 

plotSeuratReduction(seur_sce, "pca", showLegend = TRUE)

seur_sce = runSeuratFindMarkers(inSCE = seur_sce, allGroup = "Seurat_louvain_Resolution0.8")

# Fetch marker genes table
markerGenes = metadata(seur_sce)[["seuratMarkers"]] 

# Order by log fold change and p value
markerGenes = markerGenes[order(-markerGenes$avg_log2FC, markerGenes$p_val),]

markerGenes3 = markerGenes

head(markerGenes)
plotSeuratGenes(inSCE = seur_sce, useAssay = "seuratNormData", plotType = "ridge", features = metadata(seur_sce)[["seuratMarkers"]]$gene.id[9:12], groupVariable = "Seurat_louvain_Resolution0.8", ncol = 2, combine = TRUE)
FeaturePlot(seur_sce, features = features)


#cxtovary = Read10X(data.dir = "~/datasets/sample3")
#cxt_seurat = CreateSeuratObject(counts = cxtovary)

#colnames(cxt_seurat)

top100_3 = markerGenes3[1:100,1]
top100_3_1000F = marker_3_1000F[1:10,1]
top100_1 = marker_1[1:100,1]
top100_2 = marker_2[1:10,1]

count = 0
for(i in top100_2){
  print(i)
  if(i %in% top100_1){
    count = count + 1
  }
}
count