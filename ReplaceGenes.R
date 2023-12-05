library(Seurat)
library(dplyr)
library(singleCellTK)
library(ggplot2)
library(Seurat)
require("Seurat.utils")
#library(SeuratDisk)
#library(SeuratData)

library(R.utils)

cxtovary = Read10X(data.dir = "~/datasets/samples_12")
cxt_seurat = CreateSeuratObject(counts = cxtovary)
ctx_seurat = NormalizeData(object = cxt_seurat)
cxt_seurat = FindVariableFeatures(object = cxt_seurat)
cxt_seurat = ScaleData(cxt_seurat, features = rownames(cxt_seurat)) #Remove features = to increase speed but reduce number of scaled genes
cxt_seurat = RunPCA(cxt_seurat)
cxt_seurat = FindNeighbors(cxt_seurat)
cxt_seurat = FindClusters(object = cxt_seurat)
cxt_seurat = RunTSNE(cxt_seurat)

cxt_non_proc = as.SingleCellExperiment(cxt_seurat)
sce = runCellQC(cxt_non_proc, sample = NULL,
                algorithms = c("QCMetrics", "scDblFinder", "decontX"),
                seed = 12345,mitoPrefix = 'nbis')#,mitoID=mtGenes,mitoIDType = "symbol")



umap_sce = runQuickUMAP(sce, reducedDimName = "QC_UMAP",seed = 2023, sample = NULL)

#plotUMAP(umap_sce,reducedDimName = "QC_UMAP")#,colorBy = "cluster")

#FeaturePlot(umap_sce, features = "QC_UMAP")
#umap_sce = runQuickUMAP(sce, reducedDimName = "QC_UMAP", seed = 2023)
#Not working yet
#Error in .manageCellVar(inSCE, sample) : 
#Specified variable 'sample' not found in colData(inSCE)

plotRunPerCellQCResults(sce)
plotScDblFinderResults(umap_sce, reducedDimName = "QC_UMAP")
#Needs umap above to function

plotDecontXResults(umap_sce, reducedDimName = "QC_UMAP")


#reportCellQC(sce)

#hist(sce$nCount_RNA)
#hist(sce$nFeature_RNA)

saveRDS(sce, file = "sample_12_postqc.rds")