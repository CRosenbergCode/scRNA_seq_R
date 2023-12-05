#temp_seur_scaled  = ScaleData(temp_seur)#, features = rownames(temp_seur)) #Remove features = to increase speed but reduce number of scaled genes

#saveRDS(temp_seur_scaled, file ="samplepooledforheatmap.rds")

#DoHeatmap(temp_seur_scaled, features = top10$gene) + NoLegend()#features = top10$gene) + NoLegend()
library(Seurat)


temp_seur = readRDS("pooled_prescale.rds")
temp_seur_scaled_regressed_count = ScaleData(temp_seur,vars.to.regress=c("nCount_RNA"))
saveRDS(temp_seur_scaled_regressed_count, file ="samplepooledREGRESSED_Countforheatmap.rds")
