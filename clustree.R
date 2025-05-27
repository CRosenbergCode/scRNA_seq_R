#Load Libraries
library(dplyr)
library(tidyverse)
library(clustree)
library(Seurat)

#Set working directory to current directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Read in previously created seurat object to save time
seurat_ob = readRDS("seurat_clustered_1_23.RDS")

#Set a range of resolutions to examine
resolution.range <- seq(from = 0.3, to = 2.1, by = 0.2)

#Find clusters using a range of resolutions
test_res <- FindClusters(object = seurat_ob, resolution = resolution.range,useReduction = "pca",algorithm = "louvain", dims = 10)
clustree(test_res,prefix = "RNA_snn_res.")

#Find clusters based on a range of different numbers of principal components
ndims.range <- seq(from = 5, to = 50, by = 5)
test_dims=c()
for(i in ndims.range){
  test_ob = FindNeighbors(seurat_ob, dims = 1:i)
  append(test_dims,Seurat::FindClusters(object = test_ob, resolution = 0.8,useReduction = "pca",algorithm = "louvain", dims = ndims.range))
}
#Plot clusters calculated above
test_dims=FindClusters(object = seurat_ob, resolution = 0.8,useReduction = "pca",algorithm = "louvain", dims = ndims.range)
clustree(test_dims,prefix = "RNA_snn_res.")



