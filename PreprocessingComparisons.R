library(Seurat)
library(dplyr)
library(singleCellTK)
library(ggplot2)
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

geom_point(aes(x=sce$nCount_RNA,y=sce$nFeature_RNA))

ncol(sce)

saveRDS(sce, file = "sample1postqc.rds")

sce_cols = subsetSCECols(sce, colData = c("total > 600", 
                                          "detected > 300"))
#Change to umap_sce when workng
ncol(sce_cols)

cxtovary = Read10X(data.dir = "~/datasets/sample2")
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
sce2 = runCellQC(test_single, sample = NULL,
                algorithms = c("QCMetrics", "scDblFinder", "decontX"),
                seed = 12345)

#umap_sce = runQuickUMAP(sce, reducedDimName = "QC_UMAP", seed = 2023)
#Not working yet
#Error in .manageCellVar(inSCE, sample) : 
#Specified variable 'sample' not found in colData(inSCE)

plotRunPerCellQCResults(sce2)
#plotScDblFinderResults(umap_sce, reducedDimName = "QC_UMAP")
#Needs umap above to function

#reportCellQC(sce)

hist(sce2$nCount_RNA)
hist(sce2$nFeature_RNA)

geom_point(aes(x=sce2$nCount_RNA,y=sce2$nFeature_RNA))

ncol(sce2)

saveRDS(sce2, file = "sample2postqc.rds")

sce_cols2 = subsetSCECols(sce2, colData = c("total > 600", 
                                          "detected > 300"))
#Change to umap_sce when workng
ncol(sce_cols2)

cxtovary = Read10X(data.dir = "~/datasets/sample3")
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
sce3 = runCellQC(test_single, sample = NULL,
                algorithms = c("QCMetrics", "scDblFinder", "decontX"),
                seed = 12345)

#umap_sce = runQuickUMAP(sce, reducedDimName = "QC_UMAP", seed = 2023)
#Not working yet
#Error in .manageCellVar(inSCE, sample) : 
#Specified variable 'sample' not found in colData(inSCE)

plotRunPerCellQCResults(sce3)
#plotScDblFinderResults(umap_sce, reducedDimName = "QC_UMAP")
#Needs umap above to function

#reportCellQC(sce)

hist(sce3$nCount_RNA)
hist(sce3$nFeature_RNA)

type(sce3$nCount_RNA)
temp_count = as.vector(sce3$nCount_RNA)
temp_features = as.vector(sce3$nFeature_RNA)
temp_count[10]
test3 = ggplot(x=temp_count,y=temp_features) + geom_point()
#test3 = ggplot(x=sce3$nCount_RNA,y=sce3$nFeature_RNA) + geom_point()
test3

ncol(sce3)

saveRDS(sce3, file = "sample3postqc.rds")


sce_cols3 = subsetSCECols(sce3, colData = c("total > 600", 
                                          "detected > 300"))
#Change to umap_sce when workng
ncol(sce_cols3)

fil = sce$nCount_RNA[sce$nCount_RNA<4000]
fil2 = sce2$nCount_RNA[sce2$nCount_RNA<4000]
fil3 = sce3$nCount_RNA[sce3$nCount_RNA<4000]

fet = sce$nFeature_RNA[sce$nFeature_RNA<2000]
fet2 = sce2$nFeature_RNA[sce2$nFeature_RNA<2000]
fet3 = sce3$nFeature_RNA[sce3$nFeature_RNA<2000]

c3 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c1 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
c2 <- rgb(10,255,30, max = 255, alpha = 80, names = "lt.green")


hist(fil, col=c1,breaks=60,xlim=c(0, 2000)) #Pink
hist(fil2, col=c2, add=TRUE) #Green
hist(fil3, col= c3, add=TRUE) #Blue


hist(fet, col=c1,breaks=60,xlim=c(0, 2000)) #Pink
hist(fet2, col=c2, add=TRUE) #Green
hist(fet3, col= c3, add=TRUE) #Blue


ggplot(x=sce$nCount_RNA) +
  geom_violin(scale = "width", adjust = 1, width = 0.5)

print("Sample 1")
mean(sce$nCount_RNA)
sd(sce$nCount_RNA)
mean(sce$nFeature_RNA)
sd(sce$nFeature_RNA)

print("Sample 2")
mean(sce2$nCount_RNA)
sd(sce2$nCount_RNA)
mean(sce2$nFeature_RNA)
sd(sce2$nFeature_RNA)

print("Sample 3")
mean(sce3$nCount_RNA)
sd(sce3$nCount_RNA)
mean(sce3$nFeature_RNA)
sd(sce3$nFeature_RNA)

plot(sce$nCount_RNA,sce$nFeature_RNA)
plot(sce2$nCount_RNA,sce2$nFeature_RNA)
plot(sce3$nCount_RNA,sce3$nFeature_RNA)
hist(sce3$nCount_RNA,breaks=30,xlim = c(0,2000))