BiocParallel::register(BiocParallel::SerialParam())
#BiocManager::install("slingshot")


library(slingshot)
#library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)

#https://bustools.github.io/BUS_notebooks_R/slingshot.html

#annot <- readRDS("./output/neuron10k/cell_type.rds")
#mat_filtered <- readRDS("./output/neuron10k/mat_filtered.rds")

#ind <- annot$labels %in% c("NPCs", "Neurons", "OPCs", "Oligodendrocytes", 
#                           "qNSCs", "aNSCs", "Astrocytes", "Ependymal")
#cells_use <- annot$cell.names[ind]
#mat_filtered <- mat_filtered[, cells_use]

#gns <- tr2g_ensembl(species = "Mus musculus", use_gene_name = TRUE, 
#                    ensembl_version = 97)[,c("gene", "gene_name")] %>% 
#  distinct()

#seu <- CreateSeuratObject(mat_filtered) %>% 
#  SCTransform() # normalize and scale
# Add cell type annotation to metadata
#seu <- AddMetaData(seu, setNames(annot$labels[ind], cells_use), 
#                   col.name = "cell_type")

#VlnPlot(seu, c("nCount_RNA", "nFeature_RNA"), pt.size = 0.1, ncol = 1, group.by = "cell_type")

#ggplot(seu@meta.data, aes(nCount_RNA, nFeature_RNA, color = cell_type)) +
#  geom_point(size = 0.5) +
#  scale_color_brewer(type = "qual", palette = "Set2", name = "cell type") +
#  scale_x_log10() +
#  scale_y_log10() +
#  theme_bw() +
  # Make points larger in legend
#  guides(color = guide_legend(override.aes = list(size = 3))) +
#  labs(x = "Total UMI counts", y = "Number of genes detected")

seu = readRDS("rdsFiles/seurat_clustered_1_23.RDS")

seu <- RunPCA(seu, npcs = 70, verbose = FALSE)
ElbowPlot(seu, ndims = 70)


sds <- slingshot(Embeddings(seu, "umap"), clusterLabels = seu$seurat_clusters, 
                 start.clus = 4, stretch = 0)


#' Assign a color to each cell based on some value
#' 
#' @param cell_vars Vector indicating the value of a variable associated with cells.
#' @param pal_fun Palette function that returns a vector of hex colors, whose
#' argument is the length of such a vector.
#' @param ... Extra arguments for pal_fun.
#' @return A vector of hex colors with one entry for each cell.
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors <- cell_pal(seu$cluster, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(seu$seurat_clusters, hue_pal())

plot(reducedDim(sds), col = cell_colors, pch = 16, cex = 0.5)
lines(sds, lwd = 2, type = 'lineages', col = 'black')


plot(reducedDim(sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(sds, lwd = 2, type = 'lineages', col = 'black')

plot(reducedDim(sds), col = cell_colors, pch = 16, cex = 0.5)
lines(sds, lwd = 2, col = 'black')

nc <- 3
pt <- slingPseudotime(sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(sds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(sds, lwd = 2, col = 'black', type = 'lineages')
}