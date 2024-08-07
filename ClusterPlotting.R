library(Seurat)
library(EnhancedVolcano)

install.packages("EnhancedVolcano")
version("EnhancedVolcano")
devtools::install_github("kevinblighe/EnhancedVolcano")
sessionInfo()



cluster10 = c("aub","gene4776","LRR","gene5192","vasa","gene2828")
#5192 = syntaxin1,gene13019=rabx1,2828=Fcp3c,gene4776=shrb

adipocytes = c("gene4776","rabx1")
#gene13019=rabx1,gene4776=shrb

hemocytes = c("gene14629","CecA1")
#gene14629 = Karl,

temp_seur = readRDS("rdsFiles/seurat_clustered_1_23.RDS")



hem_clust = subset(x = temp_seur, subset = seurat_clusters == c(7,11))
hem_markers = FindMarkers(hem_clust, 
                      ident.1 = 7, 
                      ident.2 = 11, 
                      verbose = FALSE)
EnhancedVolcano(hem_markers, 
                rownames(hem_markers),
                x ="avg_log2FC", 
                y ="p_val_adj")#,pCutoff = 1e-02)