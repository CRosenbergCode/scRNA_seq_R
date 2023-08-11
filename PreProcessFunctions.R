library(Seurat)
library(dplyr)
library(singleCellTK)
library(ggplot2)
library(Seurat)
require("Seurat.utils")
#library(SeuratDisk)
#library(SeuratData)

preprocess3files <- function(dir_name){
  #singleCellTK()
  #cxtovary = Read10X(data.dir = "~/datasets/sample3")
  print(paste("~/datasets/",dir_name))
  cxtovary = Read10X(data.dir = paste("~/datasets/",dir_name,sep=""))
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
  #substring(dir_path)
  #saveRDS(sce, file = "sample1postqc.rds")
  saveRDS(sce, file = paste(dir_name,"postqc.rds",sep=""))
}

preprocess3files("sample1")
preprocess3files("sample2")
preprocess3files("sample3")
preprocess3files("sample3_auto")
preprocess3files("sample3_1000E")
preprocess3files("sample3_1000F")


#features_path = os.path.join(matrix_dir, "features.tsv.gz")
#feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]

# list of gene names, e.g. 'MIR1302-2HG'

library(R.utils)

#Requires R.utils library
replaceGeneName = function(featLoc, namesLoc){
  #Remove .gz for certain parts
  cutFeatLoc = substring(featLoc,1,nchar(featLoc)-3)
  feats = read.table(gzfile(featLoc))   
  annotations = read.csv(file = namesLoc, 
                         header=TRUE,skipNul = TRUE)
  
  ans = rep(TRUE, length(annotations[,4]))
  for (i in seq(annotations[,4])){
    if (nchar(annotations[i,4]) == 0){
      ans[i] = FALSE
    } 
  }
  filtered_an = annotations[ans,]
  filtered_genes = annotations[ans,1]

  nrow(feats)
  for (i in seq(nrow(feats))){
    name = feats[i,2] #1 or 2 works
    if (name %in% filtered_genes){
      temp = filtered_an[filtered_an[1] == name,4]
      if (length(temp)>1){
        temp = temp[1]
      }
      temp = gsub("\\s", "", temp)
      #print(temp)
      feats[i,1] = temp
      feats[i,2] = temp
    }
  }
  
  write.table(feats, file=cutFeatLoc, quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
  gzip(cutFeatLoc,overwrite=TRUE)
}



replaceGeneName("sample1/features.tsv.gz","Cxt_annot_corrections.csv")
replaceGeneName("sample2/features.tsv.gz","Cxt_annot_corrections.csv")
replaceGeneName("sample3/features.tsv.gz","Cxt_annot_corrections.csv")
replaceGeneName("sample3_1000E/features.tsv.gz","Cxt_annot_corrections.csv")
replaceGeneName("sample3_1000F/features.tsv.gz","Cxt_annot_corrections.csv")
replaceGeneName("sample3_auto/features.tsv.gz","Cxt_annot_corrections.csv")


sce = readRDS("sample1postqc.rds")

ncol(sce)


sce_cols = subsetSCECols(sce, colData = c("total > 75", 
                                          #                                      "detected > 20"))
                                          #sce_cols = subsetSCECols(sce, colData = c("total > 600", 
                                          "detected > 300"))
#Change to umap_sce when working
ncol(sce_cols)

seur_sce = runSeuratNormalizeData(inSCE = sce, useAssay = "decontXcounts", normAssayName = "seuratNormData", normalizationMethod = "LogNormalize", scaleFactor = 10000)

seur_sce = runSeuratFindHVG(inSCE = seur_sce, useAssay = "decontXcounts", method = "vst", hvgNumber = 2000, createFeatureSubset = "hvf")

# Print names of top 10 variable features
print(getTopHVG(inSCE = seur_sce, method = "vst", hvgNumber = 10))

plotSeuratHVG(seur_sce, labelPoints = 10)

seur_sce = runSeuratPCA(inSCE = seur_sce, useAssay = "seuratNormData", reducedDimName = "pca", nPCs = 50, seed = 42, scale = TRUE, useFeatureSubset = "hvf")

seur_sce = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 10) 

seur_sce = runSeuratFindMarkers(inSCE = seur_sce, allGroup = "Seurat_louvain_Resolution0.8")

markerGenes = metadata(seur_sce)[["seuratMarkers"]] 

markerGenes = markerGenes[order(-markerGenes$avg_log2FC, markerGenes$p_val),]

head(markerGenes)

getMarkerGenes = function(RDSfile){
  sce = readRDS(RDSfile)
  
  sce_cols = subsetSCECols(sce, colData = c("total > 75", 
                                            "detected > 20"))
  #sce_cols = subsetSCECols(sce, colData = c("total > 600", 
                                            #"detected > 300"))
  #Change to umap_sce when working

  seur_sce = runSeuratNormalizeData(inSCE = sce, useAssay = "decontXcounts", normAssayName = "seuratNormData", normalizationMethod = "LogNormalize", scaleFactor = 10000)
  
  seur_sce = runSeuratFindHVG(inSCE = seur_sce, useAssay = "decontXcounts", method = "vst", hvgNumber = 2000, createFeatureSubset = "hvf")
  
  #FindAllMarkers()
  
  seur_sce = runSeuratPCA(inSCE = seur_sce, useAssay = "seuratNormData", reducedDimName = "pca", nPCs = 50, seed = 42, scale = TRUE, useFeatureSubset = "hvf")
  
  seur_sce = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 10) 
  
  seur_sce = runSeuratFindMarkers(inSCE = seur_sce, allGroup = "Seurat_louvain_Resolution0.8")
  
  markerGenes = metadata(seur_sce)[["seuratMarkers"]] 
  
  markerGenes = markerGenes[order(-markerGenes$avg_log2FC, markerGenes$p_val),]
  
  head(markerGenes)
  
  return(markerGenes)
}

marker_3_auto = getMarkerGenes("sample3_autopostqc.rds")
#marker_3_1000E = getMarkerGenes("sample3_1000Epostqc.rds")
marker_3_1000F = getMarkerGenes("sample3_1000Fpostqc.rds")
marker_1 = getMarkerGenes("sample1postqc.rds")
marker_2 = getMarkerGenes("sample2postqc.rds")


getFullAnalysis = function(RDSFile){
  
}

compareTopGenes = function(markFiles,n){
  topVec = vector(length=length(markFiles))
  
}