library(Seurat)
library(dplyr)
library(singleCellTK)
library(dittoSeq)
library(ggplot2)
library(ggraph)
library(igraph)
library(tidygraph)
library(clustree)
library(reticulate)

#Reads data in format of the output of the 10X Genomics CellRanger Pipeline
#The directory should have at least three files: a features.tsv, a 
cxtovary = Read10X(data.dir = "~/datasets/samplePooledMito")

#The single cell experiment library is used to perform QC

cxt_seurat = CreateSeuratObject(counts = cxtovary)
#We will run some quick initial analysis and clustering for QC purposes, but will perform more careful analysis after QC 
#that will be used for analysis
#A few metrics (most especially scDoubletFinder) require initial clustering
cxt_seurat = NormalizeData(object = cxt_seurat)
cxt_seurat = FindVariableFeatures(object = cxt_seurat)
cxt_seurat = ScaleData(cxt_seurat)#, features = rownames(cxt_seurat)) #Remove features = to increase speed but reduce number of scaled genes
cxt_seurat = RunPCA(cxt_seurat)
cxt_seurat = FindNeighbors(cxt_seurat)
cxt_seurat = FindClusters(object = cxt_seurat)
cxt_seurat = RunTSNE(cxt_seurat)

#Check a single gene to see if expression is present
seur_sce['PIWIL2']

#Get the number of cells before QC
ncol(cxt_seurat)

cxt_non_proc = as.SingleCellExperiment(cxt_seurat)
#Set seed to ensure reproducible results
set.seed(12345)
#Run QC with general merics, doublet finder, and contamination
#Genes with the "MT-" prefix will be denoted as mitochondrial
sce = runCellQC(cxt_non_proc, sample = NULL,
                 algorithms = c("QCMetrics", "scDblFinder", "decontX"), mitoGeneLocation = NULL, mitoPrefix="MT-",
                 seed = 12345)

#Examples of how to look at quantiles of variables, which can be very informative for setting thresholds
#median(sce[["mito_percent"]])
#quantile(sce[["mito_percent"]],probs=seq(0,1,.1))
#quantile(sce[["mito_percent"]],probs=seq(0,1,0.01))
#quantile(sce[["nCount_RNA"]],probs=seq(0,1,.01))
#quantile(sce[["nFeature_RNA"]],probs=seq(0,1,.01))
#quantile(sce[["decontX_contamination"]],probs=seq(0,1,.01))

#Example of plotting the relationship between two variables
#ggplot(x=sce[["mito_percent"]], y=sce[["mito_detected"]]) + geom_point(x=sce[["mito_percent"]], y=sce[["mito_detected"]])
#plot(x=sce[["mito_percent"]], y=sce[["mito_detected"]])
#sce[["mito_detected"]]

#Create linear model using number of mitochondria as the independent and mito 
model_mito <- lm(sce[["mito_percent"]] ~ sce[["mito_detected"]])
#Extract coefficients from model
mito_tolerance = model_mito$coefficients[2] + sqrt(diag(vcov(model_mito)))[2]
mito_int = model_mito$coefficients[1] + sqrt(diag(vcov(model_mito)))[1]

umap_sce = runQuickUMAP(sce, reducedDimName = "QC_UMAP",seed = 2023, sample = NULL)


#Keep only cells that adhere to the following conditions:
#Called as a single cell rather than 2 or more cells by scDoubletFinder
#At least 30% of reads map to the known genome, this often varies greatly between cell types and can be as low as 55% even in humans

#Cells must have a minimum of 300 different genes expressed and 525 detected mRNA, these were based on general recommended cutoffs
#and specific empirical examination of data.

#Cells must not a mitochondrial mRNA percentage greater than their predicted percentage based on their number of mitochondria from a linear model plus twice the standard deviation of the model

#No maximum mRNA threshold was set as there appear to be biologically relevant ovarian cell types, such as nurse and stretch cells, which have extremely high mRNA expression
#For a single cell. 
sce_cols = subsetSCECols(sce, colData = c("total > 525", 
                                          "detected > 300",paste("mito_percent < ",mito_int," + ",2*mito_tolerance,"*mito_detected",sep="")
                                          ,'scDblFinder_doublet_call == "Singlet"',"decontX_contamination < 0.7"))

#Get the number of cells remaining after QC
ncol(sce_cols)

#Adjust and normalize using standard logarithmic algorithm
seur_sce = runSeuratNormalizeData(inSCE = sce_cols, useAssay = "decontXcounts", normAssayName = "seuratNormData", normalizationMethod = "LogNormalize", scaleFactor = 10000)

#Find the 5000 most variable genes, these will be used for subsequent assays such as PCA
seur_sce = runSeuratFindHVG(inSCE = seur_sce, useAssay = "decontXcounts", method = "vst", hvgNumber = 5000, createFeatureSubset = "hvf")

# Print names of top 10 variable features
#print(getTopHVG(inSCE = seur_sce, method = "vst", hvgNumber = 10))

#Plot highly variable genes and label the 10 most variable
#plotSeuratHVG(seur_sce, labelPoints = 10)

#Determine the first 50 principal components based on the above variable features, these will be used for subsequent analysis such as clustering
#This line is used for visualization purposes, we will perform this again later for subsequent analysis
seur_sce = runSeuratPCA(inSCE = seur_sce, useAssay = "seuratNormData", reducedDimName = "pca", nPCs = 50, seed = 42, scale = TRUE, useFeatureSubset = "hvf")

#Create the UMAP, which is used purely for 2D visualization to attempt to show 
#UMAP tends to preserve the relationship BETWEEN clusters while TSNE tends to better preserve relationships WITHIN clusters
#Both are used in single cell publications, I tend to prefer UMAP for our questions of interest
seur_sce = runSeuratUMAP(inSCE = seur_sce,reducedDimName = "umap",seed = 42)
#Plot the UMAP
DimPlot(temp_seur, reduction = "umap")

#A plot of the variance explained by each principal component
#The "elbow" helps us determine how many PCs to use for subsequent analysis
plotSeuratElbow(inSCE = seur_sce)

#Calculate the jackstraw plot, another way to look at PCs
seur_sce <- runSeuratJackStraw(inSCE = seur_sce, useAssay = "seuratNormData", dims = 50)
#Plot JackStraw
plotSeuratJackStraw(inSCE = seur_sce, dims = 50)


#Find clusters using a k-nearest neighbors algorithm, based on the first 10 principal components. 
seur_sce = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 10) 



#Plot first two principal components
#PCAPlot(temp_seur) + xlab("Principle Component 1") + ylab("Principle Component 2") + labs(color = "Cluster")
#Calculate and plot the TSNE, an alternative to UMAP as mentioned above
#RunTSNE(temp_seur)
#TSNEPlot(temp_seur)

#Find differentially expressed genes for each cluster
#This will only look for genes that are overexpressed 
#To look for underexpressed genes set onlyPos=FALSE
seur_sce = runSeuratFindMarkers(inSCE = seur_sce, allGroup = "Seurat_louvain_Resolution0.8",
                                minPCT = 0.25,threshUse = 0.25,onlyPos = TRUE)

#Fetch marker genes table and save to a variable
markerGenes = metadata(seur_sce)[["seuratMarkers"]] 

#Order by log fold change and p value to increase readability
markerGenes = markerGenes[order(-markerGenes$avg_log2FC, markerGenes$p_val),]

#Make ridgeplot of a specific gene, such as vitellogenin receptor
RidgePlot(object = seurat_ob,features="VGR")
#Plot a number of genes based on position in the marker gene list
plotSeuratGenes(inSCE = seur_sce, plotType = "ridge", features = metadata(seur_sce)[["seuratMarkers"]]$gene.id[17:20], groupVariable = "Seurat_louvain_Resolution0.8")#, ncol = 2, combine = TRUE)
plotSeuratGenes(inSCE = seur_sce, useAssay = "seuratNormData", plotType = "ridge", features = c("OSK"), groupVariable = "Seurat_louvain_Resolution0.8", combine = TRUE)

#Take data from a singlecellexperiment object, used for QC, to a Seurat object, used for analysis.
temp_seur = as.Seurat(seur_sce, data = "logcounts")#counts = "counts")#, data = "logcounts")

#Find clusters using a k-nearest neighbors algorithm, based on the first 10 principal components. 
temp_seur = FindNeighbors(temp_seur, reduction = "pca",dims = 1:10)
temp_seur = FindClusters(temp_seur,resolution = 0.8, algorithm = "louvain", dims = 1:10) 

#Get top 10 marker genes for each cluster
markerGenes %>%
  group_by(cluster1) %>%
  top_n(n = 10, wt = avg_log2FC) %>% 
  arrange(by_group = cluster1) -> top10OG_3


#Find differentially expressed genes for each cluster
#This will only look for genes that are overexpressed 
temp_seur.markers <- FindAllMarkers(temp_seur, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.01)

#To look for underexpressed genes set onlyPos=FALSE
temp_seur.markers_neg = FindAllMarkers(temp_seur, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.1)


#Get top 10 marker genes for each cluster

temp_seur.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

temp_seur.markers_neg %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10_neg

#Examples of subclustering, finding clusters inside of cluster
#Not utilized for current manuscript as it tended to be uninformative for mixed population clusters in our data
#sub_clusts <- FindSubCluster(temp_seur, "2", "test", subcluster.name = "unknown",  resolution = 0.75, algorithm = 1)

#DimPlot(sub_clusts, reduction = "umap", group.by = "unknown", label = TRUE, label.size = 6)

#sub_clusts = SetIdent(sub_clusts, value = sub_clusts@meta.data$unknown)

#sub_clusts.markers = FindAllMarkers(sub_clusts, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.1)

#sub_clusts.markers %>%
#  group_by(cluster) %>%
#  top_n(n = 10, wt = avg_log2FC) -> top10_sub_clust

#Plot percentage and fold change of marker genes 
#Can help determine which clusters are most easily separable based on presence or absence of a gene
ggplot(as.data.frame(top10)) + geom_point(aes(x=pct.1, y=pct.2,color=cluster)) + xlab("Percentage Expressing Inside of Cluster") + ylab("Percentage Expressing Outside of Cluster")
#
ggplot(as.data.frame(top10)) + geom_point(aes(x=pct.2, y=avg_log2FC,color=cluster)) + xlab("Percentage Expressing Outside of Cluster") + ylab("Log-fold difference in Expression")


#Compare two specific clusters, such as 10 and 11, to find marker genes between only those clusters
nurseclust = subset(x = temp_seur, subset = seurat_clusters == c(10,11) )
maybenurse.markers <- FindAllMarkers(nurseclust, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)#,ident.1=10,ident.2=12)

mean(GetAssayData(object = clust, slot = 'data')["VGR",])

head(AverageExpression(object = temp_seur, group.by = c('ident', 'groups'))$RNA)

#Obtain metadata for each cluster such as number of cells, number of unique genes expressed, and total number of mrna 
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
}

saveRDS(temp_seur, file = "pooled_prescale.rds")
#temp_seur = readRDS("sample1scaled.rds")

temp_seur=readRDS("rdsFiles/seurat_clustered_1_23.RDS")

#Prepare data for plotting heatmap of marker gene expression by cluster per cell
#Please note that this will take a long time if scaling with all genes, see below
#This process is also incredibly memory intensive if including all genes and can only be performed on a server such as cctsi
temp_seur_scaled  = ScaleData(temp_seur)#, features = rownames(temp_seur)) #Remove features = to increase speed but reduce number of scaled genes

#The same but regressing against the number of RNA in a cell, making it closter to a proportion than an absolute count
#temp_seur_scaled_regressed = ScaleData(temp_seur,vars.to.regress=c("nCount_RNA","nFeature_RNA"))
# save histogram in pdf format in current directory
#pdf(file="heatmap_sample_1_new.pdf")

#Save above data to avoid recreating it.
saveRDS(temp_seur_scaled, file ="sampleScaledForHeatmap.RDS")

#histogram we want to save
#This is a heatmap (colors indicating z-score) of the top ten marker genes of each cluster across all cells
DoHeatmap(temp_seur_scaled, features = top10$gene) #+ NoLegend()#features = top10$gene) + NoLegend()

#Dotplots of specific genes of interest by cluster
cd_genes <- c("Vgr", "Vg","CecA1")

cd_genes <- c("shrb", "spict","rabx1","aub","crq","mob2","syx1A","hml","gene5192","gene13019")

cd_genes <- c("aub","gene5192","gene13019")
labeller_cars <- c("aub"="aubergine","gene5192"="syntaxin",
                   "gene13019"="rabx1")
DotPlot(object = temp_seur, features = cd_genes) + ylab("Cluster") + xlab("Gene")+scale_x_discrete(labels=function(x) str_replace_all(x, labeller_cars))
