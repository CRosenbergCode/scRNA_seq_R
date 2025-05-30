library(Seurat)
cellClust=readRDS("rdsFiles/seurat_clustered_1_23.RDS")

#s_genes=c("BUB1","CDK1","NHP2","AURKA","NCAPD2","TOP2B","CTCF")

# Acquire the G2M phase genes        
#g2m_genes=c("CDC45","MCM5","MCM6","CHAF1B","TYMS")

#List of putative s-phase markers
s_genes=c("BUB1","CDK1","NHP2","AURKA","NCAPD2","TOP2B","CTCF",
          "gene832","gene6977","gene11612","gene10281","gene4612",
          "gene4356","gene6705","gene8169")

#List of putative G2M phase markers
g2m_genes=c("CDC45","MCM5","MCM6","CHAF1B","TYMS","gene8484")

#Obtain markers for each cell
seurat_phase <- CellCycleScoring(cellClust,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes)
#Obtain highly variable features for cell cycle
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                    selection.method = "vst",
                                    nfeatures = 501, 
                                    verbose = FALSE)
#Extract Features
hvgList = seurat_phase@assays$RNA@var.features 

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA and color by cell cycle phase
seurat_phase <- RunPCA(seurat_phase)

#Run UMAP for visualization purpases
#testUMAP = runQuickUMAP(seurat_phase, reducedDimName = "QC_UMAP",seed = 2023, sample = NULL)
testUMAP = RunUMAP(seurat_phase, dims = 1:10)

#Print metadata
for(i in seq(0,13)){
  print(i)
  tempclust=subset(x=seurat_phase,subset=seurat_clusters==i)
  print(paste("Undecided:",sum(tempclust$Phase=="Undecided")))
  print(paste("S Phase:",sum(tempclust$Phase=="S")))
  print(paste("G1 Phase:",sum(tempclust$Phase=="G1")))
  print(paste("G2/M Phase:",sum(tempclust$Phase=="G2M")))
  print(median(tempclust$total))#[["total"]]))
  print(median(tempclust$detected))#[["detected"]]))
  print(median(tempclust$total)/median(tempclust$detected))#[["total"]]))
  print(median(tempclust$mito_percent))
  print(median(tempclust$mito_detected))
  print(ncol(tempclust))
}

#Preallocate our vectors
numCells=rep(0,14)
medGenes=rep(0,14)
medUMIs=rep(0,14)
medCopies=rep(0,14)
medMitoNum=rep(0,14)
medMitoPer=rep(0,14)
medDoublet=rep(0,14)
avDoublet=rep(0,14)
cycleUndecided=rep(0,14)
cycleS=rep(0,14)
cycleG1=rep(0,14)
cycleG2M=rep(0,14)


#Calculate metadata for each cluster
for(i in seq(0,13)){
  print(i)
  j=i+1
  tempclust=subset(x=seurat_phase,subset=seurat_clusters==i)
  
  cycleUndecided[j]=sum(tempclust$Phase=="Undecided")
  cycleS[j]=sum(tempclust$Phase=="S")
  cycleG1[j]=sum(tempclust$Phase=="G1")
  cycleG2M[j]=sum(tempclust$Phase=="G2M")
  
  numCells[j]=sum(seurat_phase$seurat_clusters==i)
  #print(median(tempclust$total))
  #print(length(test[test==i]))
  medUMIs[j]=median(tempclust$total)
  medGenes[j]=median(tempclust$detected)
  medCopies[j]=median(tempclust$total)/median(tempclust$detected)
  medMitoPer[j]=median(tempclust$mito_percent)
  medMitoNum[j]=median(tempclust$mito_detected)
  medDoublet[j]=median(tempclust$scDblFinder_doublet_score)
  avDoublet[j]=mean(tempclust$scDblFinder_doublet_score)
}

clusterNum=seq(0,13)

#Create columns for dataframe to save metadata in
clustSumDF=data.frame(clusterNum,medGenes,medUMIs,medCopies,medMitoPer,medMitoNum,cycleUndecided,cycleS,cycleG1,cycleG2M,
                      medDoublet,avDoublet)
#Manually add column names
colnames(clustSumDF) = c("Cluster Number","Median # of Genes", "Median # of UMIs","Average UMIs per Gene","Median Mito. %",
                        "Median # of Mito", "# of CC Undecided","# of CC S","# of CC G1","Number of CC G2/M",
                        "Median Doublet","Average Doublet")
#Save cell cluster metadata for use in supplemental figure
write.csv(clustSumDF,file='CellClusterCharacteristics.csv',na='')

# Visualize the PCA, grouping by cell cycle phase
DimPlot(seurat_phase,
        reduction = "umap",
        group.by= "Phase")
#Grouped by cluster for comparison
DimPlot(seurat_phase,
        reduction = "umap",
        group.by= "seurat_clusters")

#Plotting lists of genes to determine if expression present
RidgePlot(cellClust, features = c("CDK1", "AURKA","MCM5", "MCM6"), ncol = 2) #No AURKA
RidgePlot(cellClust, features = c("BUB1", "NHP2","CDC45", "CHAF1B"), ncol = 2) #No BUB1
RidgePlot(cellClust, features = c("NCAPD2", "TOP2B","CTCF", "TYMS"),log=TRUE, ncol = 2)

#Plot specific potential marker genes of interest
RidgePlot(cellClust, features = c("gene6705", "gene4356","gene8169"), ncol = 2) #No BUB1
