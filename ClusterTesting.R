#install.packages("factoextra")
library(factoextra)
#library(dplyr)

head(Embeddings(temp_seur, reduction = "pca")[, 1:5])
temp_seur_1 = runUMAP(temp_seur, reducedDimName = "QC_UMAP",seed = 2023, sample = NULL)
testdf = Embeddings(temp_seur_1, reduction = "pca")

fviz_nbclust(testdf, kmeans, method = "silhouette")

set.seed(123)
km.res = kmeans(testdf, centers=6, nstart = 100)

fviz_cluster(km.res, testdf,geom = "point")


#umap_km = runQuickUMAP(km.res, reducedDimName = "QC_UMAP",seed = 2023, sample = NULL)

plotUMAP(temp_seur_1,reducedDimName = "QC_UMAP")#,colorBy = km.res)
