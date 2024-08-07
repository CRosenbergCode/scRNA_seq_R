library(Seurat)
library(dplyr)
library(singleCellTK)
library(dittoSeq)
library(ggplot2)
library(reticulate)


#Data was analyzed on a linux server and uses Unix style directory notation
#

cxtovary = Read10X(data.dir = "~/datasets/samplePooledMito")

cxt_seurat = CreateSeuratObject(counts = cxtovary)
cxt_seurat = NormalizeData(object = cxt_seurat)
cxt_seurat = FindVariableFeatures(object = cxt_seurat)
cxt_seurat = ScaleData(cxt_seurat)#, features = rownames(cxt_seurat)) #Remove features = to increase speed but reduce number of scaled genes
cxt_seurat = RunPCA(cxt_seurat)
cxt_seurat = FindNeighbors(cxt_seurat)
cxt_seurat = FindClusters(object = cxt_seurat)
cxt_seurat = RunTSNE(cxt_seurat)

#Determine number of barcodes called as cells prior to using QC
ncol(cxt_seurat)

cxt_non_proc = as.SingleCellExperiment(cxt_seurat)
set.seed(12345)
sce = runCellQC(cxt_non_proc, sample = NULL,
                 algorithms = c("QCMetrics", "scDblFinder", "decontX"), mitoGeneLocation = NULL, mitoPrefix="MT-",
                 seed = 12345)#,geneSetList = list(mtGenes),geneSetListLocation = "rownames")#,mitoID=mtGenes)#,mitoPrefix = 'MT-',mitoIDType = "symbol")#,mitoID=mtGenes,mitoIDType = "symbol")

model_mito <- lm(sce[["mito_percent"]] ~ sce[["mito_detected"]])
mito_tolerance = model_mito$coefficients[2] + sqrt(diag(vcov(model_mito)))[2]
mito_int = model_mito$coefficients[1] + sqrt(diag(vcov(model_mito)))[1]


model_rna <- lm(sce[["nCount_RNA"]] ~ sce[["nFeature_RNA"]])# + sce[["nFeature_RNA"]]:sce[["nCount_RNA"]])
rna_tolerance = model_rna$coefficients[2] + sqrt(diag(vcov(model_rna)))[2]
rna_int = model_rna$coefficients[1] + sqrt(diag(vcov(model_rna)))[1]

umap_sce = runQuickUMAP(sce, reducedDimName = "QC_UMAP",seed = 2023, sample = NULL)


#The following command filters out columns based on the following values
#

sce_cols = subsetSCECols(sce, colData = c("total > 525", 
                                          "detected > 300",paste("mito_percent < ",mito_int," + ",2*mito_tolerance,"*mito_detected",sep="")
                                          ,'scDblFinder_doublet_call == "Singlet"',"decontX_contamination < 0.7"))

ncol(sce_cols)

seur_sce = runSeuratNormalizeData(inSCE = sce_cols, useAssay = "decontXcounts", normAssayName = "seuratNormData", normalizationMethod = "LogNormalize", scaleFactor = 10000)

seur_sce = runSeuratFindHVG(inSCE = seur_sce, useAssay = "decontXcounts", method = "vst", hvgNumber = 5000, createFeatureSubset = "hvf")

seur_sce = runSeuratPCA(inSCE = seur_sce, useAssay = "seuratNormData", reducedDimName = "pca", nPCs = 50, seed = 42, scale = TRUE, useFeatureSubset = "hvf")

seur_sce = runSeuratUMAP(inSCE = seur_sce,reducedDimName = "umap",seed = 42)

plotSeuratElbow(inSCE = seur_sce)

seur_sce = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 10) 

seur_sce = runSeuratFindMarkers(inSCE = seur_sce, allGroup = "Seurat_louvain_Resolution0.8",
                                minPCT = 0.25,threshUse = 0.25,onlyPos = TRUE)


markerGenes = metadata(seur_sce)[["seuratMarkers"]] 
#
# # Order by log fold change and p value
markerGenes = markerGenes[order(-markerGenes$avg_log2FC, markerGenes$p_val),]

RidgePlot(object = seurat_ob,features="VGR")
plotSeuratGenes(inSCE = seur_sce, plotType = "ridge", features = metadata(seur_sce)[["seuratMarkers"]]$gene.id[17:20], groupVariable = "Seurat_louvain_Resolution0.8")#, ncol = 2, combine = TRUE)
plotSeuratGenes(inSCE = seur_sce, useAssay = "seuratNormData", plotType = "ridge", features = c("OSK"), groupVariable = "Seurat_louvain_Resolution0.8", combine = TRUE)

temp_seur = as.Seurat(seur_sce, data = "logcounts")#counts = "counts")#, data = "logcounts")

temp_seur = FindNeighbors(temp_seur, reduction = "pca",dims = 1:10)
temp_seur = FindClusters(temp_seur,resolution = 0.8, algorithm = "louvain", dims = 1:10) 

#saveRDS(temp_seur, file = "sample3clustered.rds")

temp_seur.markers <- FindAllMarkers(temp_seur, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)

temp_seur.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

ggplot(as.data.frame(top10)) + geom_point(aes(x=pct.1, y=pct.2,color=cluster)) + xlab("Percentage Expressing Inside of Cluster") + ylab("Percentage Expressing Outside of Cluster")

ggplot(as.data.frame(top10)) + geom_point(aes(x=pct.2, y=avg_log2FC,color=cluster)) + xlab("Percentage Expressing Outside of Cluster") + ylab("Log-fold difference in Expression")


top10['pct.1']

test = temp_seur[["seurat_clusters"]]

clust_1 = temp_seur[,test == 1]
ncol(clust_1)

reclust = 




saveRDS(temp_seur, file = "pooled_prescale.rds")

temp_seur=readRDS("rdsFiles/seurat_clustered_1_23.RDS")

temp_seur_scaled  = ScaleData(temp_seur)#, features = rownames(temp_seur)) #Remove features = to increase speed but reduce number of scaled genes

saveRDS(temp_seur_scaled, file ="sampleScaledForHeatmap.RDS")


#saveRDS(temp_seur_scaled_regressed, file ="samplepooledREGRESSEDforheatmap.rds")


# a histogram we want to save
DoHeatmap(temp_seur_scaled, features = top10$gene) #+ NoLegend()#features = top10$gene) + NoLegend()

data("pbmc_small")
cd_genes <- c("Vgr", "Vg","CecA1")

cd_genes <- c("shrb", "spict","rabx1","aub","crq","mob2","syx1A","hml","gene5192","gene13019")

cd_genes <- c("aub","gene5192","gene13019")
labeller_cars <- c("aub"="aubergine","gene5192"="syntaxin",
                   "gene13019"="rabx1")
DotPlot(object = temp_seur, features = cd_genes) + ylab("Cluster") + xlab("Gene")+scale_x_discrete(labels=function(x) str_replace_all(x, labeller_cars))


cd_genes <- c("AUB1","AUB2","AUB3","AUB4","AUB5","AUB6")
#labeller_cars <- c("aub"="aubergine","gene5192"="syntaxin",
                   #"gene13019"="rabx1")
DotPlot(object = temp_seur, features = cd_genes) + ylab("Cluster") + xlab("Gene")#+scale_x_discrete(labels=function(x) str_replace_all(x, labeller_cars))

cd_genes <- c("SGPL1","IDI1","TUBA8","MGME1")
labeller_cars <- c("aub"="aubergine","gene5192"="syntaxin",
                   "gene13019"="rabx1")
DotPlot(object = temp_seur, features = cd_genes)# + ylab("Cluster") + xlab("Gene")+scale_x_discrete(labels=function(x) str_replace_all(x, labeller_cars))
