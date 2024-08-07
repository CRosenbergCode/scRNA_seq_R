#Place these two files in same directory
seurat_ob = readRDS("rdsFiles/seurat_clustered_1_23.RDS")
#old_ob.markers = readRDS("seurat_markers_1_23.RDS")
seurat_ob.markers = readRDS("rdsFiles/seurat_markers_1_23.RDS")

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

getClusterSummary = function(clusts = seq(0,13)){
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
      exp = sum(GetAssayData(object = clust, slot = 'data')[gene,] > 0)
      #print(len)
      #print(exp)
      print(paste("Percentage of Cells: ",exp/len,sep=""))
      print(paste("Med. Expression: ",median(GetAssayData(object = clust, slot = 'data')[gene,]),sep=""))
      print(paste("Mean Expression: ",mean(GetAssayData(object = clust, slot = 'data')[gene,]),sep=""))
    }
  }
}

#makeClusterCSV = function(genes,clusts = seq(0,13)){
  retFrame = data.frame(matrix("NA", nrow = length(genes), ncol = 11))
  labs = c("Gene",">=80",">=60",">=50",">=40",">=30",">=20",">=10",">=2","<2","No Expression")
  colnames(retFrame) = labs
  retFrame[1] = genes
  count = 1
  for(gene in genes){
    addString = rep("", 10)
    for(i in clusts){
      clust = subset(x = seurat_ob, subset = seurat_clusters == i)
      exp = sum(GetAssayData(object = clust, slot = 'data')[gene,] > 0)
      clen = ncol(GetAssayData(object = clust, slot = 'data'))
      cellPercent = exp/clen
      print(i)
      print(cellPercent)
      if(exp == 0){
        addString[1] = paste(addString[1],i,sep=",")
      }
      else if(cellPercent < 0.02){
        addString[2] = paste(addString[2],i,sep=",")
      }
      else if(cellPercent < 0.1){
        addString[3] = paste(addString[3],i,sep=",")
      }
      else if(cellPercent >= 0.8){
        addString[10] = paste(addString[10],i,sep=",")
      }
      else if(cellPercent >= 0.6){
        addString[9] = paste(addString[9],i,sep=",")
      }
      else{
        addString[floor(cellPercent*10)+3] = paste(addString[floor(cellPercent)+1],i,sep=",") #First to 1 #Second to 3
      }
    }
    addString = rev(addString)

    for(i in seq(1,10)){
      if(addString[i] != ''){
        retFrame[count,i+1] = sub('^.', '', addString[i])
      }
    }
    count = count+1
  }
  return(retFrame)
} 
#Just more complicated and not worth it

makeClusterCSV = function(genes,clusts = seq(0,13)){
  retFrame = data.frame(matrix("NA", nrow = length(genes), ncol = 11))
  labs = c("Gene",">=80",">=60",">=50",">=40",">=30",">=20",">=10",">=2","<2","No Expression")
  colnames(retFrame) = labs
  retFrame[1] = genes
  count = 1
  for(gene in genes){
    addString = rep("", 10)
    for(i in clusts){
      clust = subset(x = seurat_ob, subset = seurat_clusters == i)
      exp = sum(GetAssayData(object = clust, slot = 'data')[gene,] > 0)
      clen = ncol(GetAssayData(object = clust, slot = 'data'))
      cellPercent = exp/clen
      #print(i)
      #print(cellPercent)
      if(exp == 0){
        addString[1] = paste(addString[1],i,sep=",")
      }
      else if(cellPercent >= 0.8){
        addString[10] = paste(addString[10],i,sep=",")
      }
      else if(cellPercent >= 0.6){
        addString[9] = paste(addString[9],i,sep=",")
      }
      else if(cellPercent >= 0.5){
        addString[8] = paste(addString[8],i,sep=",")
      }
      else if(cellPercent >= 0.4){
        addString[7] = paste(addString[7],i,sep=",")
      }
      else if(cellPercent >= 0.3){
        addString[6] = paste(addString[6],i,sep=",")
      }
      else if(cellPercent >= 0.2){
        addString[5] = paste(addString[5],i,sep=",")
      }
      else if(cellPercent >= 0.0001){
        addString[4] = paste(addString[4],i,sep=",")
      }
      else if(cellPercent > 0.02){
        addString[3] = paste(addString[3],i,sep=",")
      }
      else if(cellPercent > 0){
        addString[2] = paste(addString[2],i,sep=",")
      }
      else{
        print(paste("Here..",i))
      }
    }
    addString = rev(addString)
    
    for(i in seq(1,10)){
      if(addString[i] != ''){
        retFrame[count,i+1] = sub('^.', '', addString[i])
      }
    }
    count = count+1
  }
  return(retFrame)
}


clusterPers = makeClusterCSV(c("VGR","CecA1"))

clusterPers = makeClusterCSV(c("VGR","AUB2","AUB3","AUB4"))

write.csv(clusterPers,file='GenePercents.csv',na='')


#write.csv(ob_name,file='example_file.csv',na='')

#getGeneSummary()
#Desired genes can be placed inside of c() with each name surrounded by quotes and separated
#by commas
#Ex: getGeneSummary(c("getGeneSummary(c("gene1","gene2"))"))
getGeneSummary(c("FK506"))

#getDiffExpresion()
#First parameter is the minimum percentage of cells expressing a marker gene in the cluster
#Second parameter is the maximum percentage of cells outside of the given cluster expressing
#the marker gene
eighty_fifteen = getDiffExpression(0.8,0.2)
#To save results to a csv file
#First parameter is the object to save (eighty_fifteen), second is the 
#desired file name surround by quotes ('example_file.csv')
write.csv(eighty_fifteen,file='example_file.csv',na='')

#getClusterSummary(): Print summary statistics about each cluster
getClusterSummary()

eighty_fifteen %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10Diff

seurat_ob.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

heme10 = top10[top10$cluster==c(8,12,13),]
hemetest= subset(top10, cluster == c(8,12,13))

filter(cluster %in% c(8, 12, 13))

seurat_ob.markers %>%
  group_by(cluster) %>%
  filter(cluster %in% c(8, 12, 13)) %>%
  top_n(n = 30, wt = avg_log2FC) -> top10heme

tail(names(sort(table(top10heme$gene))), 5)

sort(table(top10heme['gene']),decreasing=TRUE)

clusterPers = makeClusterCSV(c("VGR"))

seurat_ob.markers %>%
  group_by(cluster) %>%
  filter(cluster %in% c(8)) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10eight

clusterPers = makeClusterCSV(top10eight$gene)
write.csv(clusterPers,file='clustereight.csv',na='')

seurat_ob.markers %>%
  group_by(cluster) %>%
  filter(cluster %in% c(10)) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10ten
clusterPers = makeClusterCSV(top10ten$gene)
write.csv(clusterPers,file='clusterten.csv',na='')

seurat_ob.markers %>%
  group_by(cluster) %>%
  filter(cluster %in% c(11)) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10eleven

clusterPers = makeClusterCSV(top10eleven$gene)
write.csv(clusterPers,file='clustereleven.csv',na='')

clustEleven = read.csv("clustereleven.csv")

seurat_ob.markers %>%
  group_by(cluster) %>%
  filter(cluster %in% c(12)) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10twelve

clusterPers = makeClusterCSV(top10twelve$gene)
write.csv(clusterPers,file='clustertwelve.csv',na='')


clusterPers = makeClusterCSV(c("CecA1","gene11937","gene10689","gene13448"))
write.csv(clusterPers,file='hemocyte.csv',na='')


old_ob.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> old10