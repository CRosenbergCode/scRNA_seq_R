library(Seurat)
library(dplyr)
library(singleCellTK)
library(dittoSeq)
library(ggplot2)

cxtovary = Read10X(data.dir = "~/datasets/samplePooledMito")
cxt_seurat = CreateSeuratObject(counts = cxtovary)
ctx_seurat = NormalizeData(object = cxt_seurat)
cxt_seurat = FindVariableFeatures(object = cxt_seurat)
cxt_seurat = ScaleData(cxt_seurat)#, features = rownames(cxt_seurat)) #Remove features = to increase speed but reduce number of scaled genes
cxt_seurat = RunPCA(cxt_seurat)
cxt_seurat = FindNeighbors(cxt_seurat)
cxt_seurat = FindClusters(object = cxt_seurat)
cxt_seurat = RunTSNE(cxt_seurat)

cxt_non_proc = as.SingleCellExperiment(cxt_seurat)
sce = runCellQC(cxt_non_proc, sample = NULL,
                algorithms = c("QCMetrics", "scDblFinder", "decontX"), mitoGeneLocation = NULL, mitoPrefix="MT-",
                seed = 12345)#,geneSetList = list(mtGenes),geneSetListLocation = "rownames")#,mitoID=mtGenes)#,mitoPrefix = 'MT-',mitoIDType = "symbol")#,mitoID=mtGenes,mitoIDType = "symbol")

plot(x=sce[["mito_percent"]], y=sce[["mito_detected"]])
sce[["mito_detected"]]

model_mito <- lm(sce[["mito_percent"]] ~ sce[["mito_detected"]])
summary(model_mito)

mito_tolerance = model_mito$coefficients[2] + sqrt(diag(vcov(model_mito)))[2]
mito_int = model_mito$coefficients[1] + sqrt(diag(vcov(model_mito)))[1]


model_rna <- lm(sce[["nCount_RNA"]] ~ sce[["nFeature_RNA"]])# + sce[["nFeature_RNA"]]:sce[["nCount_RNA"]])
summary(model_rna)

rna_tolerance = model_rna$coefficients[2] + sqrt(diag(vcov(model_rna)))[2]
rna_int = model_rna$coefficients[1] + sqrt(diag(vcov(model_rna)))[1]

sce_cols = subsetSCECols(sce, colData = c("total > 525", 
                                          "detected > 300",paste("mito_percent < ",mito_int," + ",mito_tolerance,"/mito_detected",sep="")
                                          ,'scDblFinder_doublet_call == "Singlet"',"decontX_contamination < 0.7"))

saveRDS(sce, file = "sample_pooled_postqc.rds")
sce = readRDS("sample_pooled_postqc.rds")

seur_sce = runSeuratNormalizeData(inSCE = sce, useAssay = "decontXcounts", normAssayName = "seuratNormData", normalizationMethod = "LogNormalize", scaleFactor = 10000)

seur_sce = runSeuratFindHVG(inSCE = seur_sce, useAssay = "decontXcounts", method = "vst", hvgNumber = 2000, createFeatureSubset = "hvf")

seur_sce = runSeuratPCA(inSCE = seur_sce, useAssay = "seuratNormData", reducedDimName = "pca", nPCs = 50, seed = 42, scale = TRUE, useFeatureSubset = "hvf")

seur_sce = runSeuratUMAP(inSCE = seur_sce,reducedDimName = "umap",seed = 42)

seur_sce = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 10) 

seur_sce = runSeuratFindMarkers(inSCE = seur_sce, allGroup = "Seurat_louvain_Resolution0.8")
# 
# # Fetch marker genes table
markerGenes = metadata(seur_sce)[["seuratMarkers"]] 
# 
# # Order by log fold change and p value
markerGenes = markerGenes[order(-markerGenes$avg_log2FC, markerGenes$p_val),]
#
temp_seur = as.Seurat(seur_sce, counts = "counts", data = "logcounts")

temp_seur <- FindNeighbors(temp_seur, dims = 1:10)
temp_seur <- FindClusters(temp_seur)

temp_seur.markers <- FindAllMarkers(temp_seur, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

temp_seur.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

saveRDS(cxt_seurat, file = "pooled_prescale.rds")
#temp_seur = readRDS("sample1scaled.rds")

temp_seur_scaled  = ScaleData(temp_seur)#, features = rownames(temp_seur)) #Remove features = to increase speed but reduce number of scaled genes

saveRDS(temp_seur_scaled, file ="samplepooledforheatmap.rds")

# a histogram we want to save
DoHeatmap(temp_seur_scaled, features = top10$gene) + NoLegend()#features = top10$gene) + NoLegend()