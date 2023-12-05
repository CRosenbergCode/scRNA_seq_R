####################################################SETTINGS#########################################################
#FASTA Sequences: file name
# IMPORTANT : NO UNDERSCORES IN FILE-Name and sequence names in FASTA FILE!
FileNameFasta <- 'cov-2.fasta'

#FASTA Sequences: path   [no '/' at the end of the path name]
#WINDOWS IMPORTANT: Replace the separators in the path name \ by /
PathNameFasta <- 'D:/Work/GitHub/projects/corona-fish/data/probes'
#Output directory and file names: generated automatically
name_base <-strsplit(FileNameFasta, "[.]")[[1]][1]
FileNameOutput <- paste('Probes__',name_base,sep="")
#Uncomment next line if you want to use a specific dG37 (ex:-32)
#dG37Desiree <- -32
#Use Ensemblformat to resolve multiple transcripts possibilities (sequence name must respect ENSEMBL Id)
EnsemblFormat <- TRUE
#dG37 probes Minimum score
ScoreMin <- 0.9
#probes minimum length [default 32]
TailleSondeMax <- 32
#probes maximum length [default 26]
TailleSondeMin <- 26
#mimimum number of nucleotides between two probes (end to start positions)
DistanceMinInterSonde <- 2
#Use of Composition rules described by the PNAS article [ TRUE to use the filter; FALSE to not use this filter]
PNASfilter <- TRUE
#Numbers of rules to be used (ex to use only rules 1,2 and 4 use c(1,2,4))
PNASfilterOption <- c(1,2,4)
#Use GC percentage filter
GCfilter <- TRUE
#GC composition must be between MinGC and MaxGC for the probe to be kept
MinGC <- 0.4
MaxGC <- 0.6
#Use of RepeatMasker Filter
MaskedFilter <- FALSE
#Percentage of masked nucleotides allowed by probes for the probe to be kept  [default 0.1]
MaxMaskedPercent <- 0.1
#Minimum number of probes by transcript for the transcript to be kept in final probes list  [we recommend that you have at least 24 probes]
minProbePerTranscrit <- 0
### Path to the installed REPEAT MASKER
# IMPORTANT 1 : THERE HAS TO BE A SPACE (!!!) behind the RepeatMasker command
# IMPORTANT 2 : The option "- species" specifies the species of the input sequence (valid NCBI Taxonomy species name). Commonly used names are human, mouse, rattus. For more details see Oligostan documentation
RepeatMaskerCommand <- '/usr/local/RepeatMasker/RepeatMasker -species human '  # IMPORTANT: THERE HAS TO BE A SPACE (!!!) behind the RepeatMasker command
####################################################END OF SETTINGS##################################################


library(ade4)
library(seqinr)
library(zoo)

ThedG37Min = -36
ThedG37Max = -28
ThedG37Step = 0.5
ThedG37Seq = seq(ThedG37Min,ThedG37Max,ThedG37Step)
theFLAPXSEQ = "CCTCCTAAGTTTCGAGCTGGACTCAGTG"
theFLAPYSEQ = "TTACACTCGGACCTCGTCGACATGCATT"
theFLAPZSEQ = "CCAGCTTCTAGCATCCATGCCCTATAAG"
RSEfilter <- FALSE

setwd(PathNameFasta)
if (!file.exists(FileNameOutput)){ dir.create(FileNameOutput) }
setwd(FileNameOutput)

# Save settings
#fileConn<-file("OLigostan_Settings.txt")
#writeLines(c("Oligostan smiFISH probe design","\t" ,Sys.Date( )), fileConn)
#close(fileConn)


# Create file names
paste(FileNameOutput,"_FILT_summary",".txt",collapse="",sep="") -> ResultFileName
paste(FileNameOutput,"_ALL_summary",".txt",collapse="",sep="")  -> RawProbesFileName

paste(FileNameOutput,"_FILT",".fasta",collapse="",sep="") -> ResultFastaSummaryFileName
paste(FileNameOutput,"_ALL",".fasta",collapse="",sep="")  -> RawProbesFastaSummaryFileName
paste(FileNameOutput,"_ALL_FLAPY",".fasta",collapse="",sep="")  -> RawProbesFlapYFastaSummaryFileName


getGandTInfosFromFastaReads <- function(FastaFile){
  
  unlist(strsplit(getName(read.fasta(FastaFile)),split="[|]")) -> fastafiletmp
  unlist(strsplit(fastafiletmp,split="_")) -> fastafiletmp
  length(unique(fastafiletmp[grep("ENSG",fastafiletmp)])) -> nog
  length(unique(fastafiletmp[grep("ENST",fastafiletmp)])) -> not
  
  return(c(genes=nog,transcrits=not))
}

setwd(PathNameFasta)
if (EnsemblFormat==F ){
  read.fasta(FileNameFasta,seqonly=T) -> multifastatmp
  multifasta <- list()
  for (i in 1:length(multifastatmp)){
    thetmpname <- paste("ENSG",i,"|ENST",i,sep="")
    #Reverse complement all seq to work directly from the probes point of view
    c(multifasta,list(as.SeqFastadna(rev(comp(s2c(multifastatmp[[i]]))),name=thetmpname,Annot=thetmpname)))->multifasta
  }
  rm(multifastatmp)
  c(genes=length(multifasta),transcrits=length(multifasta))->gantinfo
}

if (EnsemblFormat==T ) {
  read.fasta(FileNameFasta) -> multifastatmp
  multifasta <- list()
  for (i in 1:length(multifastatmp)){
    #Reverse complement all seq to work directly from the probes point of view
    c(multifasta,list(as.SeqFastadna(rev(comp(multifastatmp[[i]])),name=names(multifastatmp)[i],Annot=getAnnot(multifastatmp)[[i]][1])))->multifasta
  }
  rm(multifastatmp)
  getGandTInfosFromFastaReads(FileNameFasta) -> gantinfo
}
setwd(FileNameOutput)

ProbesTmp <- NULL
ProbesTmpNames <- NULL

if(!exists("dG37Desiree")){
  ThedG37SeqTmp <- ThedG37Seq
}
if(exists("dG37Desiree")){
  ThedG37SeqTmp <- dG37Desiree
}

WhichMax <- function(x){
  max(x) -> themax
  which(x==themax) -> thesize
  res <- c(thesize,themax)
  if(length(thesize)>=2) res <- c(0,themax)
  return(res)
}

#This function has been modified to use an linear score calculation and
#not an relative one.
#dG37ScoreCalc <- function(ThedG37,DesireddG = -33){
#  ThedG37 <- abs(ThedG37)
#  DesireddG <- abs(DesireddG)
#  ThedG37 - DesireddG -> tmpcalc
#  ((ThedG37[tmpcalc>=0] - DesireddG) / DesireddG) -> ThedG37[tmpcalc>=0]
#  ((DesireddG - (ThedG37[tmpcalc<0] )) / DesireddG) -> ThedG37[tmpcalc<0]
#  return(1-ThedG37)
#}

dG37ScoreCalc <- function(ThedG37,DesireddG = -33){
  valtemp <- abs(ThedG37-DesireddG)
  valnorm <- (-0.1 * valtemp) + 1
  return(valnorm)
}

ConvertRNASeq2DeltaGat37 <- function(RNASeq){
  NbBase = nchar(RNASeq)
  RNASeq <- toupper(RNASeq)
  DimVal  <- c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")
  dGVal37 <- c(-0.2,-1.5,-0.9,-1.0,-1.0,-2.2,-1.2,-1.4,-0.8,-2.4,-1.5,-1.0,-0.3,-1.4,-1.0,-0.4)
  RNADimThermoTable <- data.frame(DimVal,dGVal37)
  rep(0,NbBase-1) -> dG
  rep(RNASeq,NbBase-1) -> dim
  substr(dim,start=1:(NbBase-1),stop=2:NbBase)->dim
  for(i in 1:length(RNADimThermoTable[,1])){
    dG[which(dim==RNADimThermoTable[i,1])] <- RNADimThermoTable[i,2]
  }
  theconv <- data.frame(dim,dG)
  return(theconv)
}

dGCalc.RNA.37 <- function(RNASeq,ProbeLength=31,doiplot=F,doihist=F,SaltConc=0.115){
  AlldG <- NULL
  RNASeqConv <- ConvertRNASeq2DeltaGat37(RNASeq=RNASeq)
  AlldG <- rollapply(RNASeqConv$dG,ProbeLength-1,sum)
  AlldG <- AlldG - ((log(SaltConc)*-0.175)-.2)
  ProbeLengthseq <- paste(length(dir()),ProbeLength,sep="_")
  return(AlldG)
}

getProbesFromRNAdG37 <- function(Seq,MinSizeProbe=30,MaxSizeProbe=32,Desireddg=-33,MinScoreValue=0.9,IncBetwProb=4,doiplot=F,doihist=F){
  
  if(length(Seq)>1) paste(toupper(Seq),collapse="") -> Seq
  DiffSize = MaxSizeProbe - MinSizeProbe
  TheTmsTmp <- NULL
  dGCalc.RNA.37(Seq,ProbeLength=MaxSizeProbe,doiplot=doiplot,doihist=doihist) -> TheTmsTmp
  
  if (identical(TheTmsTmp, numeric(0))){
    return(NULL)
  }
  
  names(TheTmsTmp) <- MinSizeProbe
  NbofProbes <- length(TheTmsTmp)
  if(DiffSize > 0) {
    for(i in seq(DiffSize-1,0,-1)){
      cbind(dGCalc.RNA.37(Seq,ProbeLength=MinSizeProbe+i,doiplot=doiplot,doihist=doihist)[1:NbofProbes],TheTmsTmp) -> TheTmsTmp
    }
    colnames(TheTmsTmp) <- seq(MinSizeProbe,MaxSizeProbe,1)
  }
  dG37ScoreCalc(TheTmsTmp,Desireddg) -> TmScores
  if(DiffSize > 0) {t(apply(TmScores,1,WhichMax)) -> BestScores
    BestScores[,1]+(MinSizeProbe-1) -> BestScores[,1]}
  else BestScores <- cbind(rep(MinSizeProbe,times=length(TmScores)),TmScores)
  cbind(BestScores,seq(1:length(BestScores[,1]))) -> BestScores
  colnames(BestScores) <- c("ProbeSize","dGScore","Pos")
  rownames(BestScores) <- NULL
  BestScores[BestScores[,2]>=MinScoreValue,] -> ValidedScores
  TheProbes <- NULL
  if(length(ValidedScores)>3){
    ValidedScores[order(ValidedScores[,3]),] -> ValidedScores
    Pointeur <- 0
    
    while(Pointeur < (nchar(Seq)) ){
      ValidedScores[ValidedScores[,3]>=Pointeur,] -> ValiTmp
      if(length(ValiTmp)>=4){
        data.frame(ValiTmp,Seq=substr(Seq,start=ValiTmp[1,3],stop=(ValiTmp[1,3]+ValiTmp[1,1]-1))) -> ValiTmp
        rbind(TheProbes,ValiTmp[1,])->TheProbes
        Pointeur <- (ValiTmp[1,3]+ValiTmp[1,1]+IncBetwProb)
      }
      else {
        Pointeur=(nchar(Seq) +1)
      }
    }
  }
  return(TheProbes)
}

for(i in 1:length(multifasta)){
  for(dG37 in ThedG37SeqTmp){
    # Sequence to short
    if ( length(multifasta[[i]]) < TailleSondeMin) {next}
    
    list(getProbesFromRNAdG37(multifasta[[i]],
                              MinSizeProbe = TailleSondeMin,
                              MaxSizeProbe = TailleSondeMax,
                              Desireddg = dG37,
                              MinScoreValue = ScoreMin,
                              IncBetwProb = DistanceMinInterSonde
    )) -> ProbesTmpTmp
    
    # No probes found for sequence
    if (is.null(ProbesTmpTmp[[1]])){
      next
    }
    
    sprintf("%s_dG%.1f",getName(multifasta[[i]]),dG37) -> probesTmpTmpName
    cbind(ProbeName=rep(probesTmpTmpName,times=length(ProbesTmpTmp[[1]][,1])),ProbesTmpTmp[[1]]) -> theMatTmp
    c(ProbesTmp,ProbesTmpTmp) -> ProbesTmp
    c(ProbesTmpNames,probesTmpTmpName) -> ProbesTmpNames
    names(ProbesTmp) <- ProbesTmpNames
  }
}
rm(ProbesTmpTmp)
rm(probesTmpTmpName)
rm(theMatTmp)
rm(ProbesTmpNames)

isOk4PNASFilter <- function(TheSeq,filtertobeuse = c(1,2,3,4,5)){
  isITokwithacomp <- TRUE
  isITokwithastac <- TRUE
  isITokwithccomp <- TRUE
  isITokwithcstac <- TRUE
  isITokwithcspec <- TRUE
  if(1 %in% filtertobeuse) isITokwithacomp = isitok4aComp(TheSeq);
  if(2 %in% filtertobeuse) isITokwithastac = isitok4aStack(TheSeq);
  if(3 %in% filtertobeuse) isITokwithccomp = isitok4cComp(TheSeq);
  if(4 %in% filtertobeuse) isITokwithcstac = isitok4cStack(TheSeq);
  if(5 %in% filtertobeuse) isITokwithcspec = isitok4cSpecStack(TheSeq);
  return(isITokwithacomp & isITokwithastac & isITokwithccomp & isITokwithcstac & isITokwithcspec)
}

isOk4GCFilter <- function(TheSeq,minGC=0.45,maxGC=0.55){
  tolower(TheSeq) -> TheSeq
  theVerdict <- FALSE
  summary(TheSeq)$compo -> tmpcompo
  tmpcompo[names(tmpcompo)=="g"] -> gnb
  tmpcompo[names(tmpcompo)=="c"] -> cnb
  (gnb + cnb) / summary(TheSeq)$length -> gcComp
  if(gcComp <= maxGC & gcComp >= minGC) theVerdict <- TRUE
  return(theVerdict)
}

isitok4aComp <- function(theProbeSeq){
  tolower(theProbeSeq) -> theProbeSeq
  theVerdict <- FALSE
  if((summary(theProbeSeq)$compo[names(summary(theProbeSeq)$compo)=="a"] / summary(theProbeSeq)$length) < 0.28) theVerdict <- TRUE
  return(theVerdict)
}

isitok4aStack <- function(theProbeSeq){
  theVerdict <- FALSE
  probeSeq <- paste(toupper(theProbeSeq),collapse="")
  if(length(grep("AAAA",probeSeq))==0) theVerdict <- TRUE
  return(theVerdict)
}

isitok4cComp <- function(theProbeSeq){
  tolower(theProbeSeq) -> theProbeSeq
  theVerdict <- FALSE
  (summary(theProbeSeq)$compo[names(summary(theProbeSeq)$compo)=="c"] / summary(theProbeSeq)$length) -> cComp
  if(cComp < 0.28 & cComp > 0.22) theVerdict <- TRUE
  return(theVerdict)
}

isitok4cSpecStack <- function(theProbeSeq){
  tolower(theProbeSeq) -> theProbeSeq
  theVerdict <- FALSE
  matrix(rep(times=6,seq(1,6,1)),ncol=6,byrow=T) -> posrow
  matrix(rep(times=6,seq(0,5,1)),ncol=6,byrow=F) -> poscol
  matrix(theProbeSeq[posrow + poscol],ncol=6,byrow=FALSE) -> theprobestartmatrix
  apply(theprobestartmatrix,1,function(vectchar){
    summary(as.SeqFastadna(vectchar))$compo -> tmpcompo
    tmpcompo[names(tmpcompo)=="c"] -> tmpcnb
    return(tmpcnb / summary(as.SeqFastadna(vectchar))$length)
  }) -> thecpercent
  if(length(thecpercent[thecpercent > 0.5])==0) theVerdict <- TRUE
  return(theVerdict)
}

isitok4cStack <- function(theProbeSeq){
  theVerdict <- FALSE
  probeSeq <- paste(toupper(theProbeSeq),collapse="")
  if(length(grep("CCCC",probeSeq))==0) theVerdict <- TRUE
  return(theVerdict)
}

WriteProbes2FastaFileWithoutProbesNames <- function(ProbesTable,filetowrite,modetowrite){
  Sequences <- strsplit(tolower(ProbesTable$Seq),"")
  names <- as.character(seq(1:length(ProbesTable$Seq)))
  write.fasta(Sequences,names,file.out=filetowrite,open=modetowrite)
}

getInfosFromProbeList <- function(probeListNb,ProbeList=ProbesTmp,FastaSeqList=multifasta){
  strsplit(names(ProbeList)[probeListNb],"_")->SeqNameTmp
  SeqNameTmp[[1]][1] -> SeqName
  if("UTR" %in% SeqNameTmp[[1]]){
    which(match(SeqNameTmp[[1]],"UTR")==1) -> UTRNamePos
    SeqNameTmp[[1]][(UTRNamePos+1)] -> UTRPos
    paste(c(SeqNameTmp[[1]][1],SeqNameTmp[[1]][UTRNamePos],SeqNameTmp[[1]][(UTRNamePos+1)]),collapse="_") -> SeqName
  }
  if(!exists("UTRPos")){
    length( getSequence(FastaSeqList[[which(getName(FastaSeqList)==SeqName)]])) -> UTRPos
  }
  theEndPos <- NULL
  theStartPos <- NULL
  theInsideUTR <- NULL
  thedG37 <- NULL
  theGC <- NULL
  theGCFilter <- NULL
  thePNAS1 <- NULL
  thePNAS2 <- NULL
  thePNAS3 <- NULL
  thePNAS4 <- NULL
  thePNAS5 <- NULL
  thePNASFilter <- NULL
  theRSESeq <- NULL
  for(i in 1:length(ProbeList[[probeListNb]][,1])){
    length( getSequence(FastaSeqList[[which(getName(FastaSeqList)==SeqName)]])) -> seqlength
    (seqlength - ProbeList[[probeListNb]][i,3] + 1) -> EndPosTmp
    (EndPosTmp - ProbeList[[probeListNb]][i,1]) -> StartPosTmp
    if(as.integer(EndPosTmp)>as.integer(UTRPos)){
      1 -> InsideUTRTmp
    }
    if(as.integer(EndPosTmp)<=as.integer(UTRPos)){
      0 -> InsideUTRTmp
    }
    c(theInsideUTR,InsideUTRTmp) -> theInsideUTR
    c(theEndPos,EndPosTmp) -> theEndPos
    c(theStartPos,StartPosTmp) -> theStartPos
    c(thedG37,dGCalc.RNA.37(as.character(ProbeList[[probeListNb]][i,4]),ProbeLength=(EndPosTmp-StartPosTmp))) -> thedG37
    c(theGC,GC(s2c(as.character(ProbeList[[probeListNb]][i,4])))) -> theGC
    as.SeqFastadna(s2c(tolower(as.character(ProbeList[[probeListNb]][i,4])))) -> seqInFastaForm
    if(GCfilter){c(theGCFilter,isOk4GCFilter(seqInFastaForm,minGC=MinGC,maxGC=MaxGC)) -> theGCFilter}
    else if(!GCfilter){c(theGCFilter,1) -> theGCFilter}
    c(thePNAS1,isitok4aComp(seqInFastaForm)) -> thePNAS1
    c(thePNAS2,isitok4aStack(seqInFastaForm)) -> thePNAS2
    c(thePNAS3,isitok4cComp(seqInFastaForm)) -> thePNAS3
    c(thePNAS4,isitok4cStack(seqInFastaForm)) -> thePNAS4
    c(thePNAS5,isitok4cSpecStack(seqInFastaForm)) -> thePNAS5
    if(PNASfilter){c(thePNASFilter,isOk4PNASFilter(seqInFastaForm,filtertobeuse=PNASfilterOption)) -> thePNASFilter}
    else if(!PNASfilter){c(thePNASFilter,1) -> thePNASFilter}
    if(RSEfilter){c(theRSESeq,isOk4RSESeqFilter(seqInFastaForm,theRSESeqs=RSESeq)) -> theRSESeq}
    else if(!RSEfilter){c(theRSESeq,1) -> theRSESeq}
  }
  (thePNAS1+thePNAS2+thePNAS3+thePNAS4+thePNAS5) -> thePNASSum
  rm("UTRPos")
  return(cbind(theStartPos,theEndPos,ProbeSize=(theEndPos-theStartPos),dG37=thedG37,GCpc=theGC,GCFilter=theGCFilter,aCompFilter=thePNAS1,aStackFilter=thePNAS2,cCompFilter=thePNAS3,cStackFilter=thePNAS4,cSpecStackFilter=thePNAS5,NbOfPNAS=thePNASSum,PNASFilter=thePNASFilter,RSESeqFilter=theRSESeq,InsideUTR=theInsideUTR))
}

list()->AllProbesWithInfo
for (aprobelistnb in 1:length(ProbesTmp)){
  NbOfProbes <- length(ProbesTmp[[aprobelistnb]][,1])
  if (NbOfProbes == 0){next}
  
  strsplit(names(ProbesTmp)[aprobelistnb],"_dG") -> GAndTNamesTmp
  GAndTNamesTmp[[1]][2] -> dGBaseScore
  rep(GAndTNamesTmp[[1]][1],length(ProbesTmp[[aprobelistnb]][,1]))->ProbesNames
  rep(GAndTNamesTmp[[1]][2],length(ProbesTmp[[aprobelistnb]][,1]))->dGOpt
  
  getInfosFromProbeList(aprobelistnb,ProbeList=ProbesTmp,FastaSeqList=multifasta) -> ProbesTmpInfo
  
  # Catch case where only one probe is found, which throws an error in cbind
  if (NbOfProbes==1) {
    c(AllProbesWithInfo,list(cbind(dGOpt,ProbesNames,t(ProbesTmpInfo[,1:3]),ProbesTmp[[aprobelistnb]][,c(4,2)],t(ProbesTmpInfo[,c(4:13,15)])))) -> AllProbesWithInfo
  } else {
    c(AllProbesWithInfo,list(cbind(dGOpt,ProbesNames,ProbesTmpInfo[,1:3],ProbesTmp[[aprobelistnb]][,c(4,2)],ProbesTmpInfo[,c(4:13,15)]))) -> AllProbesWithInfo
  }
}


if(!exists("dG37Desiree")){
  AllGenePassedProbesNb <- NULL
  for (genesnb in 1:(length(AllProbesWithInfo)/length(ThedG37SeqTmp))){
    GenePassedProbes <- NULL
    for(dGnb in 1:length(ThedG37SeqTmp)){
      AllProbesWithInfo[[dGnb+(length(ThedG37SeqTmp)*(genesnb-1))]] -> AllPTmp
      length(which((AllPTmp$GCFilter + AllPTmp$PNASFilter + AllPTmp$RSESeqFilter) == 3)) -> NbOfPassedProbes
      c(GenePassedProbes,NbOfPassedProbes) -> GenePassedProbes
    }
    rbind(AllGenePassedProbesNb,GenePassedProbes) -> AllGenePassedProbesNb
  }
  (AllGenePassedProbesNb >= minProbePerTranscrit) -> AllGenePassedProbes
  apply(AllGenePassedProbes,2,sum) -> NbofTranscWithEnoughProbes
  apply(AllGenePassedProbesNb,2,sum) -> NbofProbesForAllTranscrit
  which(NbofTranscWithEnoughProbes==max(NbofTranscWithEnoughProbes)) -> theMaxRank
  which(NbofProbesForAllTranscrit==max(NbofProbesForAllTranscrit[theMaxRank])) -> theMaxRank
  if(length(theMaxRank)%%2==0){
    ThedG37Seq[theMaxRank[length(theMaxRank)/2]] -> theGooddG37Val
  }else{
    ThedG37Seq[theMaxRank[which(theMaxRank==median(theMaxRank))]] -> theGooddG37Val
  }
  
  TheSelecteddGProbes <- list()
  which(ThedG37Seq==theGooddG37Val) -> TheDesireedG37SeqNb
  for (genesnb in 1:(length(AllProbesWithInfo)/length(ThedG37Seq))){
    AllProbesWithInfo[[TheDesireedG37SeqNb+(length(ThedG37Seq)*(genesnb-1))]] -> AllPTmp
    c(TheSelecteddGProbes,list(AllPTmp)) -> TheSelecteddGProbes
  }
}

if(exists("dG37Desiree")){
  AllProbesWithInfo -> TheSelecteddGProbes
}

if(MaskedFilter){
  TmpSaveFileName <- "TheRepeatMaskerTmpFile.fasta"
  WriteProbes2FastaFileWithoutProbesNames(as.data.frame(do.call(rbind, TheSelecteddGProbes)),TmpSaveFileName,"w")
  system(paste(RepeatMaskerCommand,PathNameFasta,"/",FileNameOutput,"/",TmpSaveFileName,collapse='',sep=''),wait=T)
  if(file.exists(paste(PathNameFasta,"/",FileNameOutput,"/",TmpSaveFileName,".masked",collapse='',sep=''))){
    read.fasta(paste(PathNameFasta,"/",FileNameOutput,"/",TmpSaveFileName,".masked",collapse='',sep='')) -> maskedfastaprobes
    theNmaskedPC <- NULL
    for (indice in 1 : length(maskedfastaprobes)){
      maskedfastaprobes[[indice]] -> seqtmp
      (summary(seqtmp)$length - sum(summary(seqtmp)$compo)) / summary(seqtmp)$length -> npc
      c(theNmaskedPC,npc) -> theNmaskedPC
    }
  } else {
    as.data.frame(do.call(rbind, TheSelecteddGProbes)) -> maskedfastaprobes
    rep(0,length(maskedfastaprobes[,2])) -> theNmaskedPC
  }
  TheSelecteddGProbesWithRepeatMasker <- NULL
  aSimpleCounter <- 1
  for(TranscritNb in 1:length(TheSelecteddGProbes)){
    fromRMPC <- aSimpleCounter
    toRMPC <- (aSimpleCounter-1) + length(TheSelecteddGProbes[[TranscritNb]]$Seq)
    RepeatMaskerPC <- theNmaskedPC[fromRMPC:toRMPC]*100
    RepeatMaskerFilter <- as.integer(RepeatMaskerPC < (MaxMaskedPercent*100))
    aSimpleCounter <- toRMPC +1
    c(TheSelecteddGProbesWithRepeatMasker,list(data.frame(TheSelecteddGProbes[[TranscritNb]], RepeatMaskerPC,RepeatMaskerFilter))) -> TheSelecteddGProbesWithRepeatMasker
  }
  TheSelecteddGProbesWithRepeatMasker -> TheSelecteddGProbes
}

TheSelecteddGProbesWithSEQS <- NULL
for(TranscritNb in 1:length(TheSelecteddGProbes)){
  cbind(as.character(TheSelecteddGProbes[[TranscritNb]]$Seq),rep(theFLAPXSEQ,length(TheSelecteddGProbes[[TranscritNb]]$Seq))) -> testtmpseqX
  cbind(as.character(TheSelecteddGProbes[[TranscritNb]]$Seq),rep(theFLAPYSEQ,length(TheSelecteddGProbes[[TranscritNb]]$Seq))) -> testtmpseqY
  cbind(as.character(TheSelecteddGProbes[[TranscritNb]]$Seq),rep(theFLAPZSEQ,length(TheSelecteddGProbes[[TranscritNb]]$Seq))) -> testtmpseqZ
  HybFlpX <- apply(testtmpseqX,1,paste,collapse="")
  HybFlpY <- apply(testtmpseqY,1,paste,collapse="")
  HybFlpZ <- apply(testtmpseqZ,1,paste,collapse="")
  c(TheSelecteddGProbesWithSEQS,list(data.frame(TheSelecteddGProbes[[TranscritNb]],HybFlpX,HybFlpY,HybFlpZ))) -> TheSelecteddGProbesWithSEQS
}

TheSelecteddGProbesWithSEQS -> TheSelecteddGProbesWithSEQSTot

TheSelecteddGProbesSorted <- NULL
for(TranscritNb in 1:length(TheSelecteddGProbesWithSEQSTot)){
  c(TheSelecteddGProbesSorted,list(TheSelecteddGProbesWithSEQSTot[[TranscritNb]][
    order(-TheSelecteddGProbesWithSEQSTot[[TranscritNb]]$NbOfPNAS)
    ,])) -> TheSelecteddGProbesSorted
}
TheSelecteddGProbesSorted -> TheSelecteddGProbesWithSEQSTot

####### Write raw results
dataResultProbesWithSEQSTot <- as.data.frame(do.call(rbind, TheSelecteddGProbesWithSEQSTot))
write.table(dataResultProbesWithSEQSTot,RawProbesFileName,sep="\t",row.names=F)

# Write as FASTA - for nblast
startDum <- as.list(as.character(dataResultProbesWithSEQSTot$theStartPos))
endDum   <- as.list(as.character(dataResultProbesWithSEQSTot$theEndPos))
namesDum <- as.list(as.character(dataResultProbesWithSEQSTot$ProbesNames))
write.fasta(as.list(as.character(dataResultProbesWithSEQSTot$Seq)), paste(namesDum,startDum,endDum,sep="--"), nbchar = 60, RawProbesFastaSummaryFileName, open = "w")

# Write as FASTA with FLAP-Y - for nblast
startDum <- as.list(as.character(dataResultProbesWithSEQSTot$theStartPos))
endDum   <- as.list(as.character(dataResultProbesWithSEQSTot$theEndPos))
namesDum <- as.list(as.character(dataResultProbesWithSEQSTot$ProbesNames))
write.fasta(as.list(as.character(dataResultProbesWithSEQSTot$HybFlpY)), paste(namesDum,startDum,endDum,sep="--"), nbchar = 60, RawProbesFlapYFastaSummaryFileName, open = "w")



##### FILTER PROBES
TheFilteredProbes <- NULL
for(TranscritNb in 1:length(TheSelecteddGProbesWithSEQSTot)){
  if(MaskedFilter){
    TheSelecteddGProbesWithSEQSTot[[TranscritNb]]$GCFilter  & TheSelecteddGProbesWithSEQSTot[[TranscritNb]]$PNASFilter  & TheSelecteddGProbesWithSEQSTot[[TranscritNb]]$RepeatMaskerFilter -> AllFilter
  }
  if(!MaskedFilter){
    TheSelecteddGProbesWithSEQSTot[[TranscritNb]]$GCFilter  & TheSelecteddGProbesWithSEQSTot[[TranscritNb]]$PNASFilter -> AllFilter
  }
  length(TheSelecteddGProbesWithSEQSTot[[TranscritNb]][AllFilter,]$Seq) -> NbOfGoodProbe
  if(NbOfGoodProbe >= minProbePerTranscrit){c(TheFilteredProbes,list(TheSelecteddGProbesWithSEQSTot[[TranscritNb]][AllFilter,])) -> TheFilteredProbes}
  rm(NbOfGoodProbe)
}

# Write filtered results

if(length(TheFilteredProbes) == 0) {
  write(x="Not probes for specified sequence found after filtering. Change filtering parameters or minimum number of probes per transcript",file = ResultFileName)
  
} else {
  dataResultProbesFiltered <- as.data.frame(do.call(rbind, TheFilteredProbes))
  write.table(dataResultProbesFiltered,ResultFileName,sep="\t",row.names=F)
  
  # Write as FASTA - for nblast
  startDum <- as.list(as.character(dataResultProbesFiltered$theStartPos))
  endDum   <- as.list(as.character(dataResultProbesFiltered$theEndPos))
  namesDum <- as.list(as.character(dataResultProbesFiltered$ProbesNames))
  write.fasta(as.list(as.character(dataResultProbesFiltered$Seq)), paste(namesDum,startDum,endDum,sep="--"), nbchar = 60, ResultFastaSummaryFileName, open = "w")
  
}

cat('\n\n\n=== Oligonstan terminated succesfully. See result folder for identified probes.\n',PathNameFasta)

if(MaskedFilter)system("rm TheRepeatMaskerTmpFile.*",wait=T)
rm(list = ls())