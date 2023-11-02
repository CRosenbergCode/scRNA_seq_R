library(Seurat)
library(dplyr)
library(singleCellTK)
library(ggplot2)
library(Seurat)
require("Seurat.utils")
#library(SeuratDisk)
#library(SeuratData)

library(R.utils)

#Requires R.utils library
#Renames from generic gene names (gene####) to gene symbols
#featLoc is the relative path of the features.tsv.gz
#Ex: "sample1/features.tsv.gz"
#namesLoc is the relative path of the file containing the names of putative genes
#Ex: "Cxt_annot_corrections.csv"
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
      feats[i,1] = temp
      feats[i,2] = temp
    }
  }
  
  write.table(feats, file=cutFeatLoc, quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
  gzip(cutFeatLoc,overwrite=TRUE)
}

#Create SingleCell experiment object from three out files from cell ranger
#Saves object as RDS for further use
#dir_name is the relative path of the directory containing the 3 files
#ex: 
preprocess3files <- function(dir_name){
  #cxtovary = Read10X(data.dir = "~/datasets/sample3")
  cxtovary = Read10X(data.dir = paste("~/datasets/",dir_name,sep=""))
  cxt_seurat = CreateSeuratObject(counts = cxtovary)
  norm_seurat = NormalizeData(object = cxt_seurat)
  var_feats = FindVariableFeatures(object = norm_seurat)
  scales = ScaleData(var_feats)
  pca = RunPCA(scales)
  neighs = FindNeighbors(pca) #Reduction = PCA? 
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
  
  #substring(dir_path)
  #saveRDS(sce, file = "sample1postqc.rds")
  saveRDS(sce, file = paste(dir_name,"postqc.rds",sep=""))
}

runSeuratProcess = function(sceQC,){
  #sce_cols = subsetSCECols(sce, colData = c("total > 75", 
  #                                      "detected > 20"))
  #Change to umap_sce when working
  #ncol(sce_cols)
  
  seur_sce = runSeuratNormalizeData(inSCE = sceQC, useAssay = "decontXcounts", normAssayName = "seuratNormData", normalizationMethod = "LogNormalize", scaleFactor = 10000)
  
  seur_sce = runSeuratFindHVG(inSCE = seur_sce, useAssay = "decontXcounts", method = "vst", hvgNumber = 2000, createFeatureSubset = "hvf")
  
  # Print names of top 10 variable features
  #print(getTopHVG(inSCE = seur_sce, method = "vst", hvgNumber = 10))
  
  plotSeuratHVG(seur_sce, labelPoints = 10)
  
  seur_sce = runSeuratPCA(inSCE = seur_sce, useAssay = "seuratNormData", reducedDimName = "pca", nPCs = 50, seed = 42, scale = TRUE, useFeatureSubset = "hvf")
  seur_sce = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 10) 
  return(seur_sce)
}

#Example Usage
preprocess3files("samples_1_2")


#features_path = os.path.join(matrix_dir, "features.tsv.gz")
#feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]

# list of gene names, e.g. 'MIR1302-2HG'

replaceGeneName("pooledE3000/features.tsv.gz","Cxt_annot_corrections_10_30.csv")
replaceGeneName("pooledE4000/features.tsv.gz","Cxt_annot_corrections_10_30.csv")

#replaceGeneName("sample1/features.tsv.gz","Cxt_annot_mod.csv")
#replaceGeneName("sample2/features.tsv.gz","Cxt_annot_mod.csv")
#replaceGeneName("sample3/features.tsv.gz","Cxt_annot_mod.csv")
#replaceGeneName("sample3_1000E/features.tsv.gz","Cxt_annot_corrections.csv")
#replaceGeneName("sample3_1000F/features.tsv.gz","Cxt_annot_corrections.csv")
#replaceGeneName("sample3_auto/features.tsv.gz","Cxt_annot_corrections.csv")


sce = readRDS("sample1postqc.rds")

ncol(sce)


sce_cols = subsetSCECols(sce, colData = c("total > 75", 
                                          #                                      "detected > 20"))
                                          #sce_cols = subsetSCECols(sce, colData = c("total > 600", 
                                          "detected > 300"))
#Change to umap_sce when working
ncol(sce_cols)

#Scaling will significantly increase the time required to run the function
seurat_processing = function(scale=FALSE,){
  seur_sce = as.Seurat(seur_sce, counts = "counts", data = "logcounts")
  
  seur_sce = runSeuratNormalizeData(inSCE = sce, useAssay = "decontXcounts", normAssayName = "seuratNormData", normalizationMethod = "LogNormalize", scaleFactor = 10000)
  
  seur_sce = runSeuratFindHVG(inSCE = seur_sce, useAssay = "decontXcounts", method = "vst", hvgNumber = 2000, createFeatureSubset = "hvf")
  
  seur_sce = runSeuratPCA(inSCE = seur_sce, useAssay = "seuratNormData", reducedDimName = "pca", nPCs = 50, seed = 42, scale = TRUE, useFeatureSubset = "hvf")
  
  seur_sce = runSeuratUMAP(inSCE = seur_sce,reducedDimName = "umap",seed = 42)
  
  seur_sce = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 10) 
  
  #temp_seur <- FindNeighbors(temp_seur, dims = 1:10)
  #temp_seur <- FindClusters(temp_seur)
  
  
  #markerGenes = metadata(seur_sce)[["seuratMarkers"]] 
  
  #markerGenes = markerGenes[order(-markerGenes$avg_log2FC, markerGenes$p_val),]
  
  #head(markerGenes)
  
  temp_seur.markers <- FindAllMarkers(seur_seur, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  temp_seur.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
  
  
  if(scale){
    temp_seur_scaled  = ScaleData(temp_seur)#, features = rownames(temp_seur)) #Remove features = to increase speed but reduce number of scaled genes
  }
  
  if(verbose){
    print(head(markerGenes))
  }
  if(save){
    saveRDS(markerGenes,file = paste(sam_name,"_markergenes",sep=""))
  }
  return(c(markerGenes,top10))
}
# Print names of top 10 variable features
print(getTopHVG(inSCE = seur_sce, method = "vst", hvgNumber = 10))

plotSeuratHVG(seur_sce, labelPoints = 10)

seur_sce = runSeuratPCA(inSCE = seur_sce, useAssay = "seuratNormData", reducedDimName = "pca", nPCs = 50, seed = 42, scale = TRUE, useFeatureSubset = "hvf")

seur_sce = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 10) 

seur_sce = runSeuratFindMarkers(inSCE = seur_sce, allGroup = "Seurat_louvain_Resolution0.8")

markerGenes = metadata(seur_sce)[["seuratMarkers"]] 

markerGenes = markerGenes[order(-markerGenes$avg_log2FC, markerGenes$p_val),]

head(markerGenes)

temp_seur.markers <- FindAllMarkers(temp_seur, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

temp_seur.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10


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

#Returns 
getAllMarkers = function(RDSfile,sam_name, verbose=FALSE,topN=10){
  if(!getRDS(sam_name,"postqc")){
    seur_sce = readRDS(RDSfile)
  }
  
  seur_sce = runSeuratPCA(inSCE = seur_sce, useAssay = "seuratNormData", reducedDimName = "pca", nPCs = 50, seed = 42, scale = TRUE, useFeatureSubset = "hvf")
  
  seur_sce = runSeuratFindClusters(inSCE = seur_sce, useReduction = "pca", resolution = 0.8, algorithm = "louvain", dims = 10) 
  
  seur_sce = runSeuratFindAllMarkers(object = seur_sce)#runSeuratFindMarkers(inSCE = seur_sce, allGroup = "Seurat_louvain_Resolution0.8")
  
  #runSeuratFindMar
  markerGenes = metadata(seur_sce)[["seuratMarkers"]] 
  
  markerGenes = markerGenes[order(-markerGenes$avg_log2FC, markerGenes$p_val),]
  
  temp_seur.markers <- FindAllMarkers(temp_seur, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  temp_seur.markers %>%
    group_by(cluster) %>%
    top_n(n = topN, wt = avg_log2FC) -> topGenes
  
  
  if(verbose){
    print(head(markerGenes))
  }
  if(save){
    saveRDS(markerGenes,file = paste(sam_name,"_markergenes",sep=""))
  }
  return(c(markerGenes,topGenes))
}

getAllMarkers("sample1postqc.rds")

compareTopGenes = function(markFiles,n){
  topVec = vector(length=length(markFiles))
  
}

#Requires R.utils library
#Renames from generic gene names (gene####) to gene symbols
#featLoc is the relative path of the features.tsv.gz
#Ex: "sample1/features.tsv.gz"
#namesLoc is the relative path of the file containing the names of putative genes
#Ex: "Cxt_annot_corrections.csv"
replaceGeneNameFitz = function(featLoc, namesLoc){
  cutFeatLoc = substring(featLoc,1,nchar(featLoc)-3)
  feats = read.table(gzfile(featLoc))   
  annotations = read.csv(file = namesLoc, 
                         header=TRUE,skipNul = TRUE)
  
  ans = rep(TRUE, length(annotations[,4]))
  for (i in seq(annotations[,4])){
    temp = toString(annotations[i,4])
    if (temp == 0){
      ans[i] = FALSE
    }
    if(temp == "NA" || temp == "-"){
      ans[i] = FALSE
    }
  }
  filtered_an = annotations[ans,]
  filtered_genes = annotations[ans,2]
  
  nrow(feats)
  for (i in seq(nrow(feats))){
    name = feats[i,2] #1 or 2 works
    if (name %in% filtered_genes){
      temp = filtered_an[filtered_an[2] == name,4]
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

replaceGeneNameFitz("sample1/features.tsv.gz","geneID_eggnog_merge_cxt_CDS.csv")
replaceGeneNameFitz("sample2/features.tsv.gz","geneID_eggnog_merge_cxt_CDS.csv")
replaceGeneNameFitz("sample3/features.tsv.gz","geneID_eggnog_merge_cxt_CDS.csv")
replaceGeneName("features.tsv.gz","Cxt_annot_corrections_8_15.csv")


dynamic_cell_filter = function(sce){
  
  return(sce)
}

getHeatmap = function(sam_name){
  getRDS = function(sam_name,stage_name){
    if(file.exists(paste(sam_name,stage_name,sep=""))){
      return(readRDS(paste(sam_name,stage_name,sep="")))
    }
    return(FALSE)
  }
  if(!getRDS(sam_name,"forheatmap")){
    if(!getRDS(sam_name,"scal")){#LOOK AT SAMPLE NAME
      
    }
  }
  return(DoHeatmap(temp_seur_scaled, features = top10$gene) + NoLegend())
}

runQC = function(sce,decont=TRUE,doublet=TRUE,mito=TRUE,plots=TRUE,prints=TRUE,quantile=0.05,RDS=FALSE){
  
  if(plots){
    plotRunPerCellQCResults(sce)
    plotScDblFinderResults(umap_sce, reducedDimName = "QC_UMAP")
    #Needs umap above to function
    
    plotDecontXResults(umap_sce, reducedDimName = "QC_UMAP")
  }
  
  if(static){
    
  }
  else{
    count = 
      features = 
      doublet = TRUE
    decont = 0.9
  }
  return(subsetted)
}

filterCells = function(sce,decont=TRUE,doublet=TRUE,mito=TRUE,plots=TRUE,RDS=FALSE,static=FALSE,staticVals=){
  #sce_cols = subsetSCECols(sce, colData = c("total > 600", 
  #                                          "detected > 300",))
  #retCells = sce[doublet<0.9]
  
  if(static){
    
  }
  else{
    count = 
      features = 
      doublet = TRUE
    decont = 0.9
  }
  return(subsetted)
  
  sce = subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 & decontX_contamination == "Singlet" & scDblFinder_doublet_call < 0.9 & nCount_RNA < 200) 
}

plotGene = function(gene,all=TRUE,ridge=TRUE,umap=TRUE,){
  plotSeuratGenes(inSCE = seur_sce, useAssay = "seuratNormData", plotType = "ridge", features = gene, groupVariable = "Seurat_louvain_Resolution0.8", ncol = 2, combine = TRUE) + ylab("Cluster Identity")
  #FeaturePlot(seur_sce, features = c("CecA1","VgR"))#features)
  FeaturePlot(seur_sce, features = gene, min.cutoff = 1, max.cutoff = 3)
  VlnPlot(pbmc3k.final, features = c(gene))
  
}


