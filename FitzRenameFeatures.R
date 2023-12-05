#require("Seurat.utils")
library(R.utils)


#Requires R.utils library
#Renames from generic gene names (gene####) to gene symbols
#featLoc is the relative path of the features.tsv.gz
#Ex: "sample1/features.tsv.gz"
#namesLoc is the relative path of the file containing the names of putative genes
#Ex: "Cxt_annot_corrections.csv"
replaceGeneNameFitz = function(featLoc, namesLoc){
  cutFeatLoc = substring(featLoc,1,nchar(featLoc)-3)
  feats = read.table(gzfile(featLoc))   
  annotations = read.csv(file = namesLoc, 
                         header=TRUE,skipNul = TRUE)
  
  ans = rep(TRUE, length(annotations[,4]))
  for (i in seq(annotations[,4])){
    temp = toString(annotations[i,4])
    if (temp == 0){
      ans[i] = FALSE
    }
    if(temp == "NA" || temp == "-"){
      ans[i] = FALSE
    }
  }
  filtered_an = annotations[ans,]
  filtered_genes = annotations[ans,2]
  
  nrow(feats)
  for (i in seq(nrow(feats))){
    name = feats[i,2] #1 or 2 works
    if (name %in% filtered_genes){
      temp = filtered_an[filtered_an[2] == name,4]
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

#ExampleUsage
replaceGeneNameFitz("sample1/features.tsv.gz","geneID_eggnog_merge_cxt_CDS.csv")
