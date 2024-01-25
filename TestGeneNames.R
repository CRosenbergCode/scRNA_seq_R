
getDuplicateGenes = function(){
  
}



replaceGeneName = function(featLoc, namesLoc){
  #Remove .gz for certain parts
  cutFeatLoc = substring(featLoc,1,nchar(featLoc)-3)
  feats = read.table(gzfile(featLoc))   
  annotations = read.csv(file = namesLoc, 
                         header=TRUE,skipNul = TRUE)
  
  ans = rep(TRUE, length(annotations[,4]))
  for (i in seq(annotations[,4])){
    if (nchar(annotations[i,4]) == 0){
      ans[i] = FALSE
    } 
  }
  filtered_an = annotations[ans,]
  filtered_genes = annotations[ans,1]
  
  nrow(feats)
  for (i in seq(nrow(feats))){
    name = feats[i,2] #1 or 2 works
    if (name %in% filtered_genes){
      temp = filtered_an[filtered_an[1] == name,4]
      if (length(temp)>1){
        temp = temp[1]
      }
      temp = gsub("\\s", "", temp)
      feats[i,1] = temp
      feats[i,2] = temp
    }
  }
  
  write.table(feats, file=cutFeatLoc, quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
  gzip(cutFeatLoc,overwrite=TRUE)
}

replaceGeneName("samplePooledMito/features.tsv.gz","Cxt_annot_corrections_1_5.csv")
