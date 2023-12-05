#Place these two files in same directory
seurat_ob = readRDS("seurat_clustered_11_28.RDS")
seurat_ob.markers = readRDS("seurat_markers_11_28.RDS")

#Helper Function
extractNumeric = function(seur_ob){
  temp = c()
  for(i in seur_ob[,]){
    temp = c(temp,i)
  }
  return(temp)
}

getDiffExpression = function(min_clust,max_other){
  filt = seurat_ob.markers[seurat_ob.markers["pct.1"]>min_clust,]
  #pct.2 is the name of the column
  #max_other is the maximum threshold
  filt2 = filt[filt["pct.2"]<max_other,]
  return(filt2)
}

getClusterSummary = function(clusts = seq(0,12)){
  for(i in clusts){
    clust = subset(x = seurat_ob, subset = seurat_clusters == i)
    print(paste("Cluster ",i,sep=""))
    print(paste("Med. RNA per Cell: ",median(extractNumeric(clust[["total"]])),sep=""))
    print(paste("Med. Genes per Cell ",median(extractNumeric(clust[["detected"]])),sep=""))
    print(paste("Med.  ",median(extractNumeric(clust[["detected"]]))/median(extractNumeric(clust[["total"]])),sep=""))
    print(paste("Med. Percent Mitochondria ",median(extractNumeric(clust[["mito_percent"]])),sep=""))
    print(paste("Number of Cells: ",ncol(GetAssayData(object = clust, slot = 'data')),sep=""))
  }
}

getGeneSummary = function(genes,clusts = seq(0,13)){
  for(gene in genes){
    for(i in clusts){
      clust = subset(x = seurat_ob, subset = seurat_clusters == i)
      print(paste("Gene: ",gene,sep=""))
      print(paste("Cluster ",i,sep=""))
      len = ncol(GetAssayData(object = clust, slot = 'data'))
      exp = sum(GetAssayData(object = clust, slot = 'data')["VGR",] > 0)
      #print(len)
      #print(exp)
      print(paste("Percentage of Cells: ",exp/len,sep=""))
      print(paste("Med. Expression: ",median(GetAssayData(object = clust, slot = 'data')[gene,]),sep=""))
      print(paste("Mean Expression: ",mean(GetAssayData(object = clust, slot = 'data')[gene,]),sep=""))
    }
  }
}


#

#write.csv(ob_name,file='example_file.csv',na='')

#getGeneSummary()
#Desired genes can be placed inside of c() with each name surrounded by quotes and separated
#by commas
#Ex: getGeneSummary(c("getGeneSummary(c("gene1","gene2"))"))
getGeneSummary(c("VGR"))


#getDiffExpresion()
#First parameter is the minimum percentage of cells expressing a marker gene in the cluster
#Second parameter is the maximum percentage of cells outside of the given cluster expressing
#the marker gene
eighty_fifteen = getDiffExpression(0.8,0.15)
#To save results to a csv file
#First parameter is the object to save (eighty_fifteen), second is the 
#desired file name surround by quotes ('example_file.csv')
write.csv(eighty_fifteen,file='example_file.csv',na='')

#getClusterSummary(): Print summary statistics about each cluster
getClusterSummary()

seurat_ob.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10