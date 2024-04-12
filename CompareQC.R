library(Seurat)
library(SeuratObject)
library(dplyr)
library(singleCellTK)
library(dittoSeq)
library(ggplot2)
library(stringr)
library(irlba)
library(Matrix)

#https://stackoverflow.com/questions/25114771/glpk-no-such-file-or-directory-error-when-trying-to-install-r-package

library(reticulate)


#remove.packages('irlba')

#install.packages('irlba')

#remove.packages('irlba')

#install.packages('irlba')

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

create_sce = function(sampleDir,fullScale = TRUE,save=FALSE,plots=FALSE){
  cxtovary = Read10X(data.dir = paste("~/datasets/",sampleDir,sep=""))
  cxt_seurat = CreateSeuratObject(counts = cxtovary)
  cxt_seurat = NormalizeData(object = cxt_seurat)
  cxt_seurat = FindVariableFeatures(object = cxt_seurat)
  if(fullScale){
    cxt_seurat  = ScaleData(cxt_seurat, features = rownames(cxt_seurat)) #Remove features = to increase speed but reduce number of scaled genes
  }
  else{
    cxt_seurat  = ScaleData(cxt_seurat) #Remove features = to increase speed but reduce number of scaled genes 
  }
  cxt_seurat  = RunPCA(cxt_seurat)
  cxt_seurat  = FindNeighbors(cxt_seurat)
  cxt_seurat  = FindClusters(object = cxt_seurat)
  cxt_seurat  = RunTSNE(cxt_seurat)
  
  cxt_sce = as.SingleCellExperiment(cxt_seurat)
  if(save){
    saveRDS(sce, file = paste(sampleDir,"postqc.rds",sep=""))
  }
  return(cxt_sce)
}

sample1sce = create_sce("AlternateCellranger/Lib1E1800",fullScale=FALSE)
sample2sce = create_sce("AlternateCellranger/Lib2E1800",fullScale=FALSE)

#for (file in c("samp1matrix.m.gz", "samp2matrix.m.gz")){
for (file in c("Lib1E1800", "Lib2E1800")){
  seurat_data <- Read10X(data.dir = paste0("~/datasets/AlternateCellranger/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

Lib1E1800@meta.data

merged_seurat <- merge(x = Lib1E1800, 
                       y = Lib2E1800, 
                       add.cell.id = c("one", "two"))



merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
metadata <- merged_seurat@meta.data
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^one_"))] <- "one"
metadata$sample[which(str_detect(metadata$cells, "^two_"))] <- "two"

merged_seurat@meta.data <- metadata

saveRDS(merged_seurat, file="merged_filtered_seurat.rds")


metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~sample)

metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)


metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)


metadata_clean <- filtered_seurat@meta.data


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

#umap_merged = RunUMAP(merged_seurat, reducedDimName = "QC_UMAP",nn.name = "weighted.nn")
merged_seurat@neighbors
umap_sce = runQuickUMAP(sce, reducedDimName = "QC_UMAP",seed = 2023, sample = NULL)


seur_sce['PIWIL2']

ncol(cxt_seurat)

#cxt_single = as.SingleCellExperiment(tsnes)
#sample1PreQC = as.SingleCellExperiment(sample1sce)
#sample2PreQC = as.SingleCellExperiment(sample2sce)
mergedSCE = as.SingleCellExperiment(merged_seurat)
set.seed(12345)
sample1QC = runCellQC(sample1sce, sample = NULL,
                algorithms = c("QCMetrics", "scDblFinder", "decontX"), mitoGeneLocation = NULL, mitoPrefix="MT-",
                seed = 12345)#,geneSetList = list(mtGenes),geneSetListLocation = "rownames")#,mitoID=mtGenes)#,mitoPrefix = 'MT-',mitoIDType = "symbol")#,mitoID=mtGenes,mitoIDType = "symbol")

sample2QC = runCellQC(sample2sce, sample = NULL,
                      algorithms = c("QCMetrics", "scDblFinder", "decontX"), mitoGeneLocation = NULL, mitoPrefix="MT-",
                      seed = 12345)#,geneSetList = list(mtGenes),geneSetListLocation = "rownames")#,mitoID=mtGenes)#,mitoPrefix = 'MT-',mitoIDType = "symbol")#,mitoID=mtGenes,mitoIDType = "symbol")

#FeaturePlot(pbmc3k.final, features = features)



sample1QCU = runQuickUMAP(sample1QC, reducedDimName = "QC_UMAP",seed = 2023, sample = NULL)

plotRunPerCellQCResults(sample1QC)
plotScDblFinderResults(sample1QCU, reducedDimName = "QC_UMAP")
plotDecontXResults(sample1QCU, reducedDimName = "QC_UMAP")


sample2QCU = runQuickUMAP(sample2QC, reducedDimName = "QC_UMAP",seed = 2023, sample = NULL)

plotRunPerCellQCResults(sample2QC)
plotScDblFinderResults(sample2QCU, reducedDimName = "QC_UMAP")
plotDecontXResults(sample2QCU, reducedDimName = "QC_UMAP")



#Needs umap above to function



#reportCellQC(sce)

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


model_rna <- lm(sce[["nCount_RNA"]] ~ sce[["nFeature_RNA"]])# + sce[["nFeature_RNA"]]:sce[["nCount_RNA"]])
#summary(model_rna)

rna_tolerance = model_rna$coefficients[2] + sqrt(diag(vcov(model_rna)))[2]
rna_int = model_rna$coefficients[1] + sqrt(diag(vcov(model_rna)))[1]


#merged.srt[["percent_mt"]] <- PercentageFeatureSet(merged.srt, pattern = "^MT-")


#testsce = subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 & decontX_contamination == "Singlet" & scDblFinder_doublet_call < 0.9 & nCount_RNA < 200)

umap_sce = runQuickUMAP(sample1QC, reducedDimName = "QC_UMAP",seed = 2023, sample = NULL)

umap_sce = runQuickUMAP(sample2QC, reducedDimName = "QC_UMAP",seed = 2023, sample = NULL)

colData(sample1QC)

plotUMAP(sample1QC,reducedDimName = "QC_UMAP",runUMAP = TRUE)#,colorBy = "metadata$sample")
plotUMAP(sample1QC,reducedDimName = "QC_UMAP",colorBy = "cluster")



split_seurat <- SplitObject(merged_seurat, split.by = "sample")

split_seurat <- split_seurat[c("one", "two")]

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
}

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
seurat_integrated = RunPCA(seurat_integrated,reduction.name = "pca2")
#runSeuratPCA(inSCE = seur_sce, useAssay = "seuratNormData", reducedDimName = "pca", nPCs = 50, seed = 42, scale = TRUE, useFeatureSubset = "hvf")
#seurat_integrated@assays
test_int = runSeuratPCA(inSCE = seurat_integrated, reducedDimName = "pca", nPCs = 50, seed = 42, scale = TRUE, useFeatureSubset = "hvf")

seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca2")

# Plot UMAP                             
DimPlot(seurat_integrated)
PCAPlot(seurat_integrated,
        split.by = "sample")  
DimPlot(seurat_integrated,
        split.by = "sample") 

#paste(sam_name,stage_name,sep="")

#testsce = subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5 & decontX_contamination == "Singlet" & scDblFinder_doublet_call < 0.9 & nCount_RNA < 200)




sce_cols = subsetSCECols(sce, colData = c("total > 525", 
                                          "detected > 300",paste("mito_percent < ",mito_int," + ",2*mito_tolerance,"*mito_detected",sep="")
                                          ,'scDblFinder_doublet_call == "Singlet"',"decontX_contamination < 0.7"))
