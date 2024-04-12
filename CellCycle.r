library(Seurat)
cellClust=readRDS("rdsFiles/seurat_clustered_1_23.RDS")

s_genes=c("BUB1","CDK1","NHP2","AURKA","NCAPD2","TOP2B","CTCF")

# Acquire the G2M phase genes        
g2m_genes=c("CDC45","MCM5","MCM6","CHAF1B","TYMS")


seurat_phase <- CellCycleScoring(cellClust,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes)

seurat_phase <- FindVariableFeatures(seurat_phase, 
                                    selection.method = "vst",
                                    nfeatures = 501, 
                                    verbose = FALSE)

hvgList = seurat_phase@assays$RNA@var.features 

zz <- file("hvg.txt", "wb")
writeBin( paste(hvgList, collapse="\n"), zz ) 
close(zz)


# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA and color by cell cycle phase
seurat_phase <- RunPCA(seurat_phase)

#testUMAP = runQuickUMAP(seurat_phase, reducedDimName = "QC_UMAP",seed = 2023, sample = NULL)
testUMAP = RunUMAP(seurat_phase, dims = 1:10)

sum(seurat_phase$Phase=="Undecided")

# Visualize the PCA, grouping by cell cycle phase
DimPlot(seurat_phase,
        reduction = "umap",
        group.by= "Phase")

RidgePlot(cellClust, features = c("CDK1", "AURKA","MCM5", "MCM6"), ncol = 2) #No AURKA
RidgePlot(cellClust, features = c("BUB1", "NHP2","CDC45", "CHAF1B"), ncol = 2) #No BUB1
RidgePlot(cellClust, features = c("NCAPD2", "TOP2B","CTCF", "TYMS"),log=TRUE, ncol = 2)

#R value for these genes?
