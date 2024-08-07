sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves,nknots=15)

saveRDS(sce,file="GAMfitSCE.RDS")