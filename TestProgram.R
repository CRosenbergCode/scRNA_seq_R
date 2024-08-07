library(Seurat)
library(dplyr)
library(singleCellTK)
library(dittoSeq)
library(ggplot2)

require(remotes)
install_version("igraph", version = "1.5.1", repos = "http://cran.us.r-project.org")

#install.packages("clustree")
install.packages("ggraph")
#install.packages("tidygraph")
remove.packages('igraph')
install.packages("igraph")
install.packages("Rglpk")

library(ggraph)
library(igraph)
library(tidygraph)
library(clustree)
library()

#https://stackoverflow.com/questions/25114771/glpk-no-such-file-or-directory-error-when-trying-to-install-r-package
#export LD_LIBRARY_PATH=/home/camcorey/GLPK/lib and export CPATH=/home/camcorey/GLPK/include

library(reticulate)
#use_python("/home/camcorey/.local/lib/python3.8/")
#Check for your local system
#For leidenalg used in clustering, updated version of 
#Must install through python with "pip install leidenalg"
#library(SeuratDisk)
#library(SeuratData)
#packageDescription("singleCellTK")$Version
#packageDescription("SummarizedExperiment")$Version
R.version
#singleCellTK()
#cxtovary = Read10X(data.dir = "~/datasets/pooledE4000")
cxtovary = Read10X(data.dir = "~/datasets/samplePooledMito")
#cxtovary = Read10X(data.dir = "~/datasets/oldSamp")

#Features fil 31 instead of 53, 22 short)

cxt_seurat = CreateSeuratObject(counts = cxtovary)
cxt_seurat = NormalizeData(object = cxt_seurat)
cxt_seurat = FindVariableFeatures(object = cxt_seurat)
cxt_seurat = ScaleData(cxt_seurat)#, features = rownames(cxt_seurat)) #Remove features = to increase speed but reduce number of scaled genes
cxt_seurat = RunPCA(cxt_seurat)
cxt_seurat = FindNeighbors(cxt_seurat)
cxt_seurat = FindClusters(object = cxt_seurat)
cxt_seurat = RunTSNE(cxt_seurat)
#cxt_seurat[["percent_mt"]] <- PercentageFeatureSet(cxt_seurat, pattern = "^MT-")
#median(cxt_seurat[["percent_mt"]])
#Version(cxt_seurat)
#cxt_seurat[["percent_mt"]]

#saveRDS(cxt_seurat, file = "sample1scaled.rds")
#tsnes = readRDS("sample1scaled.rds")

#mtGenes <- readLines("IntermediateCTarsalisMitoGenes.txt")
#mtGenes
#'MT-ATP6'

#order(test_cxt_seurat[["RNA"]]@meta.features[['vst.mean']], # Sequence of vectors of the same length
#      decreasing = TRUE, # Whether to sort in increasing or decreasing order
#      na.last = TRUE,)     # Whether to put NA values at the beginning or at the end
#sort.list(test_cxt_seurat[["RNA"]]@meta.features[['vst.mean']], decreasing = TRUE)

#Get Generally most highly expressed genes

#test_cxt_seurat[["RNA"]]@meta.features[sort.list(test_cxt_seurat[["RNA"]]@meta.features[['vst.mean']], decreasing = TRUE),]

#SaveH5Seurat(tsnes, overwrite = TRUE)
#DimPlot(tsnes, reduction = "tsne")
#cxt_non_proc['para']
#cxt_non_proc['nbis-gene-2']
#cxt_non_proc['Aub1']
#cxt_non_proc['MT-nbis-gene-3']
#cxt_non_proc['ND4L']


#No ATP6, ATP8,para, COX1,2,3,
#Yes ND4L, ND6,
#cxt_non_proc['AUB']
#cxt_non_proc['AUB1']
#cxt_non_proc['AUB2']
#cxt_non_proc['AUB3']
#cxt_non_proc['AUB4']
#cxt_non_proc['AUB5']
#cxt_non_proc['AUB']

#cxt_single['SUMO1']
#cxt_single['PIWIL2']

seur_sce['PIWIL2']

ncol(cxt_seurat)

#cxt_single = as.SingleCellExperiment(tsnes)
cxt_non_proc = as.SingleCellExperiment(cxt_seurat)
set.seed(12345)
sce = runCellQC(cxt_non_proc, sample = NULL,
                 algorithms = c("QCMetrics", "scDblFinder", "decontX"), mitoGeneLocation = NULL, mitoPrefix="MT-",
                 seed = 12345)#,geneSetList = list(mtGenes),geneSetListLocation = "rownames")#,mitoID=mtGenes)#,mitoPrefix = 'MT-',mitoIDType = "symbol")#,mitoID=mtGenes,mitoIDType = "symbol")

#sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
#median(sce[["mito_percent"]])
#quantile(sce[["mito_percent"]],probs=seq(0,1,.1))
#quantile(sce[["mito_percent"]],probs=seq(0,1,0.01))
#quantile(sce[["nCount_RNA"]],probs=seq(0,1,.01))
#quantile(sce[["nFeature_RNA"]],probs=seq(0,1,.01))
#quantile(sce[["decontX_contamination"]],probs=seq(0,1,.01))

#ggplot(x=sce[["mito_percent"]], y=sce[["mito_detected"]]) + geom_point(x=sce[["mito_percent"]], y=sce[["mito_detected"]])
#plot(x=sce[["mito_percent"]], y=sce[["mito_detected"]])
#sce[["mito_detected"]]

model_mito <- lm(sce[["mito_percent"]] ~ sce[["mito_detected"]])
#summary(model_mito)

#model_mito$coefficients[2]
#test = summary(model_mito)#$Std.error
#sqrt(diag(vcov(model_mito)))[2]
mito_tolerance = model_mito$coefficients[2] + sqrt(diag(vcov(model_mito)))[2]
mito_int = model_mito$coefficients[1] + sqrt(diag(vcov(model_mito)))[1]

#Test

model_rna <- lm(sce[["nCount_RNA"]] ~ sce[["nFeature_RNA"]])# + sce[["nFeature_RNA"]]:sce[["nCount_RNA"]])
#summary(model_rna)

rna_tolerance = model_rna$coefficients[2] + sqrt(diag(vcov(model_rna)))[2]
rna_int = model_rna$coefficients[1] + sqrt(diag(vcov(model_rna)))[1]


#merged.srt[["percent_mt"]] <- PercentageFeatureSet(merged.srt, pattern = "^MT-")


#testsce = subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 & decontX_contamination == "Singlet" & scDblFinder_doublet_call < 0.9 & nCount_RNA < 200)

umap_sce = runQuickUMAP(sce, reducedDimName = "QC_UMAP",seed = 2023, sample = NULL)

#plotUMAP(umap_sce,reducedDimName = "QC_UMAP")#,colorBy = "cluster")


#paste(sam_name,stage_name,sep="")



sce_cols = subsetSCECols(sce, colData = c("total > 525", 
                                          "detected > 300",paste("mito_percent < ",mito_int," + ",2*mito_tolerance,"*mito_detected",sep="")
                                          ,'scDblFinder_doublet_call == "Singlet"',"decontX_contamination < 0.7"))

#test_doublet = subsetSCECols(sce, colData = c('scDblFinder_doublet_call == "Singlet"'))
#ncol(test_doublet)

#test_mito = subsetSCECols(sce, colData = c(paste("mito_percent < ",mito_int," + ",2*mito_tolerance,"*mito_detected",sep="")))
                                        
#print(mito_int)
#print(mito_tolerance)

ncol(sce_cols)
#Was 2738 without naming, lets change

#sce_cols = subsetSCECols(sce, colData = c("total > 600", 
                                          #"detected > 300","mito_percent < 10",'scDblFinder_doublet_call == "Singlet"',"decontX_contamination < 0.9"))


#FeaturePlot(umap_sce, features = "QC_UMAP")
#umap_sce = runQuickUMAP(sce, reducedDimName = "QC_UMAP", seed = 2023)
#Not working yet
#Error in .manageCellVar(inSCE, sample) : 
#Specified variable 'sample' not found in colData(inSCE)

#library(DoubletFinder)

#plotRunPerCellQCResults(sce)
#plotScDblFinderResults(umap_sce, reducedDimName = "QC_UMAP")
#Needs umap above to function

#plotDecontXResults(umap_sce, reducedDimName = "QC_UMAP")


#reportCellQC(sce)

#hist(sce$nCount_RNA)
#hist(sce$nFeature_RNA)

#saveRDS(sce, file = "sample_pooled_postqc.rds")
#sce = readRDS("sample_pooled_postqc.rds")

#Sam_name refers the 
#For example, 
#stage_name refers to the stage of the process you are in
#For 

#getSampleName = function(){
  
#}

#getRDS = function(sam_name,stage_name){
#  if(file.exists(paste(sam_name,stage_name,sep=""))){
#    return(paste(sam_name,stage_name,sep=""))
#  }
#  return(FALSE)
#}

#sce = readRDS("sample1postqc.rds")

#ncol(sce)

#sce_cols = subsetSCECols(sce, colData = c("total > 75", 
#                                      "detected > 20"))
#sce_cols = subsetSCECols(sce, colData = c("total > 600", 
#                                          "detected > 300","mito_percent < 10",'decontX_contamination == "Singlet"',"scDblFinder_doublet_call < 0.7"))
#sce_cols = subsetSCECols(sce, colData = c("total > 600", 
#                                          "detected > 300","mito_percent < 10",'scDblFinder_doublet_call == "Singlet"',"decontX_contamination < 0.9"))


#testsce = subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & mito_percent < 5 & decontX_contamination == "Singlet" & scDblFinder_doublet_call < 0.9 & nCount_RNA < 200)

#Use cutoff of less than two standard deviations away

#

#Change to umap_sce when working
#ncol(sce_cols)

seur_sce = runSeuratNormalizeData(inSCE = sce_cols, useAssay = "decontXcounts", normAssayName = "seuratNormData", normalizationMethod = "LogNormalize", scaleFactor = 10000)

seur_sce = runSeuratFindHVG(inSCE = seur_sce, useAssay = "decontXcounts", method = "vst", hvgNumber = 5000, createFeatureSubset = "hvf")

# Print names of top 10 variable features
#print(getTopHVG(inSCE = seur_sce, method = "vst", hvgNumber = 10))

#plotSeuratHVG(seur_sce, labelPoints = 10)

seur_sce = runSeuratPCA(inSCE = seur_sce, useAssay = "seuratNormData", reducedDimName = "pca", nPCs = 50, seed = 42, scale = TRUE, useFeatureSubset = "hvf")

seur_sce = runSeuratUMAP(inSCE = seur_sce,reducedDimName = "umap",seed = 42)

plotSeuratElbow(inSCE = seur_sce)

#seur_sce <- runSeuratJackStraw(inSCE = seur_sce, useAssay = "seuratNormData", dims = 50)

#Plot JackStraw
#plotSeuratJackStraw(inSCE = seur_sce, dims = 50)

seur_sce = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 10) 
#seur_sce_5 = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 5) 
#seur_sce_15 = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 15) 

#seur_sce_2 = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 2) 
#seur_sce_8 = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 8) 
#seur_sce_12 = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 12) 

#DimPlot(in_seur_sce)


#10:16
#5:15
#15:17
#2:20
#8:16
#12:16

#seur_sce = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 10) 

#pdf(file = 'output/my_plot.pdf')
#par(mar = c(4.1, 4.4, 4.1, 1.9), xaxs="i", yaxs="i")
#plot(flowers$weight, flowers$shootarea, 
#       xlab = "weight (g)",
#       ylab = expression(paste("shoot area (cm"^"2",")")),
#       xlim = c(0, 30), ylim = c(0, 200), bty = "l",
#       las = 1, cex.axis = 0.8, tcl = -0.2,
#       pch = 16, col = "dodgerblue1", cex = 0.9)
#text(x = 28, y = 190, label = "A", cex = 2)
#dev.off()

#pdf(file = 'testOut/my_plot.pdf')

#plotSeuratReduction(seur_sce, "pca", showLegend = TRUE)
#dev.off()

#plotSeuratReduction(seur_sce, "umap", showLegend = TRUE)
#FeaturePlot(seur_sce, features = "CecA1", min.cutoff = 1, max.cutoff = 3)

#plotSeuratGenes(inSCE = seur_sce, useAssay = "seuratNormData", plotType = "ridge", features = "MT-nbis-gene-2", groupVariable = "Seurat_louvain_Resolution0.8", ncol = 2, combine = TRUE) + ylab("Cluster Identity")

#FeaturePlot(temp_seur,features = pca)

#DimPlot(temp_seur) + xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + labs(color = "Cluster")
                                                                    
#PCAPlot(temp_seur) + xlab("Principle Component 1") + ylab("Principle Component 2") + labs(color = "Cluster")

#RunTSNE(temp_seur)

#TSNEPlot(temp_seur)

#sce_2 <- findMarkerDiffExp(seur_sce, useAssay = "logcounts", method = "wilcox",
#                          cluster = "cluster",
#                          log2fcThreshold = 0, fdrThreshold = 0.05,
#                          minClustExprPerc = 0, maxCtrlExprPerc = 1,
#                          minMeanExpr = 0)
#markerHm <- plotFindMarkerHeatmap(seur_sce, topN = 5, log2fcThreshold = 0, 
#                               fdrThreshold = 0.05, minClustExprPerc = 0.5, 
#                               maxCtrlExprPerc = 0.4, minMeanExpr = 0, 
#                               rowLabel = TRUE)

#sub_clustering()

seur_sce = runSeuratFindMarkers(inSCE = seur_sce, allGroup = "Seurat_louvain_Resolution0.8",
                                minPCT = 0.25,threshUse = 0.25,onlyPos = TRUE)
# 
# # Fetch marker genes table
markerGenes = metadata(seur_sce)[["seuratMarkers"]] 
# 
# # Order by log fold change and p value
markerGenes = markerGenes[order(-markerGenes$avg_log2FC, markerGenes$p_val),]
# 
#markerGenes1 = markerGenes
# 
# head(markerGenes)
RidgePlot(object = seurat_ob,features="VGR")
plotSeuratGenes(inSCE = seur_sce, plotType = "ridge", features = metadata(seur_sce)[["seuratMarkers"]]$gene.id[17:20], groupVariable = "Seurat_louvain_Resolution0.8")#, ncol = 2, combine = TRUE)
plotSeuratGenes(inSCE = seur_sce, useAssay = "seuratNormData", plotType = "ridge", features = c("OSK"), groupVariable = "Seurat_louvain_Resolution0.8", combine = TRUE)
#plotSeuratGenes(inSCE = seur_sce, useAssay = "seuratNormData", plotType = "ridge", features = c("MT-ND4","RPS27","MDH1","RPL9"), groupVariable = "Seurat_louvain_Resolution0.8", ncol = 2, combine = TRUE)


#FeaturePlot(seur_sce, features = c("CecA1","VgR"))#features)
#FeaturePlot(temp_seur_scaled, features = "SREBF2")#features)

#cxtovary = Read10X(data.dir = "~/datasets/sample3")
#cxt_seurat = CreateSeuratObject(counts = cxtovary)

#colnames(cxt_seurat)

# top100_3 = markerGenes3[1:100,1]
# top100_3_1000F = marker_3_1000F[1:10,1]
# top100_1 = marker_1[1:100,1]
# top100_2 = marker_2[1:10,1]
# 
# count = 0
# for(i in top100_2){
#   print(i)
#   if(i %in% top100_1){
#     count = count + 1
#   }
# }
# count

temp_seur = as.Seurat(seur_sce, data = "logcounts")#counts = "counts")#, data = "logcounts")

temp_seur = FindNeighbors(temp_seur, reduction = "pca",dims = 1:10)
temp_seur = FindClusters(temp_seur,resolution = 0.8, algorithm = "louvain", dims = 1:10) 

#temp_seur = runSeuratFindClusters(temp_seur,useReduction = "pca",resolution = 0.8, algorithm = "louvain", dims = 10) 


#cluster.name =

#test_temp = FindClusters(temp_seur,resolution = 0.8, algorithm = 4) 


#temp_seur <- FindClusters(temp_seur)

#temp_seur = clusts

markerGenes %>%
  group_by(cluster1) %>%
  top_n(n = 10, wt = avg_log2FC) %>% 
  arrange(by_group = cluster1) -> top10OG_3


#saveRDS(temp_seur, file = "sample3clustered.rds")

temp_seur.markers <- FindAllMarkers(temp_seur, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.01)

#head(temp_seur.markers)

#filt = temp_seur.markers[temp_seur.markers["pct.1"]>0.85,]

#head(filt)

#filt2 = filt[filt["pct.2"]<0.15,]

#head(filt2)

#write.csv(filt2,file='differential_85_15.csv',na='')

#write.csv(top10,file='top10.csv',na='')


#Subcluster stuff

#FeaturePlot(temp_seur, features = "aub") + xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + labs(color = "Average Expression")

#DimPlot(temp_seur, reduction = "umap", label = TRUE, label.size = 6 )

#temp_seur = FindClusters(temp_seur,graph.name="test",resolution = 0.8, algorithm = "louvain", dims = 10) 

#temp_seur = runSeuratFindClusters(temp_seur,resolution = 0.8, algorithm = "louvain", dims = 10) 

#temp_seur <- FindNeighbors(temp_seur,graph.name="test", dims = 1:10)
#sub_clusts <- FindSubCluster(temp_seur, "2", "test", subcluster.name = "unknown",  resolution = 0.75, algorithm = 1)
#0.5 = 3 clusters for 2
#0.
#DimPlot(sub_clusts, reduction = "umap", group.by = "unknown", label = TRUE, label.size = 6)

#sub_clusts = SetIdent(sub_clusts, value = sub_clusts@meta.data$unknown)

#sub_clusts.markers = FindAllMarkers(sub_clusts, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.1)

#sub_clusts.markers %>%
#  group_by(cluster) %>%
#  top_n(n = 10, wt = avg_log2FC) -> top10_sub_clust

temp_seur.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

temp_seur.markers_neg = FindAllMarkers(temp_seur, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.1)

temp_seur.markers_neg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10_neg

#temp_seur.markers_len = FindAllMarkers(temp_seur, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)

#temp_seur.markers_len %>%
#  group_by(cluster) %>%
#  top_n(n = 10, wt = avg_log2FC) -> top10_len

#temp_seur.markers_per <- FindAllMarkers(temp_seur, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#temp_seur.markers_per %>%
#  group_by(cluster) %>%
#  top_n(n = 10, wt = avg_log2FC) -> top10

#write.csv(top10_len,file='top10GenesClusters11_5.csv',na='')


#CLUSTREE


clustree(temp_seur, prefix = "K")
paste0(assay, "_snn_res.")


clustree_overlay(nba_clusts, prefix = "K", x_value = "PC1", y_value = "PC2")




#Totem










ggplot(as.data.frame(top10)) + geom_point(aes(x=pct.1, y=pct.2,color=cluster)) + xlab("Percentage Expressing Inside of Cluster") + ylab("Percentage Expressing Outside of Cluster")

ggplot(as.data.frame(top10)) + geom_point(aes(x=pct.2, y=avg_log2FC,color=cluster)) + xlab("Percentage Expressing Outside of Cluster") + ylab("Log-fold difference in Expression")


top10['pct.1']

#temp_seur.markers %>%
  #group_by(cluster) %>%
   #nrows()

#nrow(temp_seur[,temp_seur[["seurat_clusters"]]==10])

test = temp_seur[["seurat_clusters"]]

#length(test[test == 10])
#nrow(temp_seur[["seurat_clusters"]])

clust_1 = temp_seur[,test == 1]
ncol(clust_1)

reclust = 

#temp_seur[["seurat_clusters"]]==10

mean_umis = 
mean_mito =
mean_feats =


  sce_cols = subsetSCECols(sce, colData = c("total > 525", 
                                            "detected > 300",paste("mito_percent < ",mito_int," + ",2*mito_tolerance,"*mito_detected",sep="")
                                            ,'scDblFinder_doublet_call == "Singlet"',"decontX_contamination < 0.7"))

clust = subset(x = temp_seur, subset = seurat_clusters == 0)

nurseclust = subset(x = temp_seur, subset = seurat_clusters == c(10,11) )


maybenurse.markers <- FindAllMarkers(nurseclust, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)#,ident.1=10,ident.2=12)


gsc1clust = subset(x = temp_seur, subset = seurat_clusters == c(10,11) )


gsc1.markers <- FindAllMarkers(gsc1clust, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)#,ident.1=10,ident.2=12)



print(clust[["RNA"]]["VgR"])
'para'

length(clust[["RNA"]]["VgR"])

testob = clust[["RNA"]]#[]
head(x = rownames(x = clust))
clust[["RNA"]][,"VGR"]

mean(GetAssayData(object = clust, slot = 'data')["VGR",])

head(AverageExpression(object = temp_seur, group.by = c('ident', 'groups'))$RNA)

saveRDS(temp_seur,file="seurat_clustered_1_23.RDS")

temp_seur=readRDS(file="rdsFiles/seurat_clustered_1_23.RDS")

saveRDS(temp_seur.markers,file="seurat_markers_1_23.RDS")

test10 = readRDS(file = "seurat_markers_11_28.RDS")

quantile(t(temp_seur[["mito_percent"]]))

ncol(GetAssayData(object = clust, slot = 'data'))
ncol(ckust)

median(as.numeric(temp_seur[["total"]][4,]))
extractNumeric = function(seur_ob,subset="total"){
  temp = c()
  for(i in temp_seur[[subset]][,]){
    temp = c(temp,i)
  }
  return(temp)
}
median(extractNumeric(clust,subset="detected"))

temp_seur[["seurat_clusters"]]

for(i in seq(0,13)){
  #clust = seur_sce[,test == i]
  clust = subset(x = temp_seur, subset = seurat_clusters == i)
  print("Test")
  print(paste("Cluster ",i,sep=""))
  #print(length(test[test==i]))
  print(median(clust[["total"]]))
  print(median(clust[["detected"]]))
  print(median(clust[["detected"]])/median(clust[["total"]]))
  print(median(clust[["mito_percent"]]))
  print(length(test[test == i]))
  #print()
}

clust[["total"]][1]

umis_per_feat = 

#temp_seur.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)

#top10


saveRDS(temp_seur, file = "pooled_prescale.rds")
#temp_seur = readRDS("sample1scaled.rds")

temp_seur=readRDS("rdsFiles/seurat_clustered_1_23.RDS")

temp_seur_scaled  = ScaleData(temp_seur)#, features = rownames(temp_seur)) #Remove features = to increase speed but reduce number of scaled genes

#temp_seur_scaled_regressed = ScaleData(temp_seur,vars.to.regress=c("nCount_RNA","nFeature_RNA"))
# save histogram in pdf format in current directory
#pdf(file="heatmap_sample_1_new.pdf")

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

Dot

cd_genes <- c("AUB1","AUB2","AUB3","AUB4","AUB5","AUB6")
#labeller_cars <- c("aub"="aubergine","gene5192"="syntaxin",
                   #"gene13019"="rabx1")
DotPlot(object = temp_seur, features = cd_genes) + ylab("Cluster") + xlab("Gene")#+scale_x_discrete(labels=function(x) str_replace_all(x, labeller_cars))

cd_genes <- c("SGPL1","IDI1","TUBA8","MGME1")
labeller_cars <- c("aub"="aubergine","gene5192"="syntaxin",
                   "gene13019"="rabx1")
DotPlot(object = temp_seur, features = cd_genes)# + ylab("Cluster") + xlab("Gene")+scale_x_discrete(labels=function(x) str_replace_all(x, labeller_cars))


#DotPlot9

#kirre
#gene6166                                                          377    0.0     6
#ID=gene14430                                                         257    0.0     7
#ID=gene1590                                                          158    1e-128  7
#ID=gene3156                                                          174    1e-111  5
#ID=gene6355                                                          141    9e-109  8
#ID=gene7631                                                          136    1e-104  9
#ID=gene8378
cd_genes <- c("gene6166","gene14430","CLINT1","DCAF10","gene6355","gene7631")

DotPlot(object = temp_seur, features = cd_genes)# + ylab("Cluster") + xlab("Gene")+scale_x_discrete(labels=function(x) str_replace_all(x, labeller_cars))


#Kwon paper genes
#SCRB3

cd_genes <- c("SCRB3","gene452","SPARC","osk","vas","gene3564")
#Believed to be oenocytoid marker]
#check SCRB9 as well
""
#Kwon paper
#SPARC universal, in 12, 13, 8, 5
DotPlot(object = temp_seur, features = cd_genes)# + ylab("Cluster") + xlab("Gene")+scale_x_discrete(labels=function(x) str_replace_all(x, labeller_cars))


#call this function to save the file 
#dev.off()

#temp_seur_scaled[["cluster"]]

#temp_seur[["seuratMarker"]] %>% group_by(cluster)

#[["seuratMarkers"]]

#temp_seur_scaled %>% group_by(Idents(temp_seur_scaled))
#markers_list = Idents(temp_seur_scaled)
#length(markers_list)
#dims(temp_seur_scaled[""])
# pdf("my_plot.pdf",         # File name
#     width = 8, height = 7, # Width and height in inches
#     bg = "white",          # Background color
#     colormodel = "cmyk"    # Color model (cmyk is required for most publications)
#     paper = "A4")          # Paper size
# 
# # Creating a plot
# plot(rnorm(20))
# 
# # Closing the graphical device
# dev.off() 

# DoHeatmap(object = temp_seur)
# 
# dittoHeatmap(seur_sce, genes,
#              annot.by = "clustering")
# 
# genes <- getGenes(seur_sce)[1:100]
# 
# plotSCEHeatmap(seur_sce,)
# #THIS IS THE GOOD ONE?
# 
sce <- runScranSNN(seur_sce, useReducedDim = "PCA", clusterName = "cluster", algorithm = "louvain", k = 4)
# 
sce <- findMarkerDiffExp(sce, useAssay = "logcounts", method = "wilcox",
                          cluster = "cluster",
                          log2fcThreshold = 0, fdrThreshold = 0.05,
                          minClustExprPerc = 0, maxCtrlExprPerc = 1,
                          minMeanExpr = 0)
# 
markerHm <- plotMarkerDiffExp(sce, topN = 5, log2fcThreshold = 0, 
                               fdrThreshold = 0.05, minClustExprPerc = 0.5, 
                               maxCtrlExprPerc = 0.4, minMeanExpr = 0, 
                               rowLabel = TRUE)
# 
# markerHm
#ggsave("test_heatmap_pic.png", plot = last_plot(), device = png(), scale = 1, width = 20, height = 20, dpi = 300)
