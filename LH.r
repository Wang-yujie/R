timestart<-Sys.time()

library(optparse)
option_list = list(
  make_option(c("-id", "--sampleId"), type="character", default=NULL, 
              help="sample ID", metavar="character"),
  make_option(c("-o", "--outputDir"), type="character", default=NULL, 
              help="location of output directory", metavar="character"),
  make_option(c("-M", "--MultFactor"), type="numeric", default=NULL,
              help="GermLineUniqMappedReads/TumorUniqMappedReads ", metavar="character"),
  make_option(c("-mc", "--minCoverageFilter"), type="numeric", default=30,
              help="minimum coverage at mismatch site [default= %default]", metavar="character"),
  make_option(c("-tpl", "--tumorPloidy"), type="numeric", default=2,
              help="tumorPloidy [default= %default]", metavar="character"),
  make_option(c("-tpu", "--tumorPurity"), type="numeric", default=0.8,
              help="tumorPurity [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt)


sampleID <- opt$sampleId
workDir  <- opt$outputDir
minCoverageFilter <- opt$minCoverageFilter
tumorPloidy <- opt$tumorPloidy
tumorPurity <- opt$tumorPurity
MultFactor <- opt$MultFactor

if (is.null(opt$sampleId) | is.null(opt$outputDir) | is.null(opt$MultFactor)){
  print_help(opt_parser)
  stop("Missing arguments.\n", call.=FALSE)  
}



require(seqinr, quietly = TRUE)
require(Biostrings, quietly = TRUE) 
require(beeswarm, quietly = TRUE)
require(zoo, quietly = TRUE) ##Y
require(Rsamtools, quietly = TRUE)

binSize<-150
gamma<-1

funCalcN_withBAF <- function(logRSites,bafSites,tumorPloidy,tumorPurity,gamma)
{
  nA <- (tumorPurity-1+bafSites*2^(logRSites/gamma)*((1-tumorPurity)*2+tumorPurity*tumorPloidy))/tumorPurity
  nB <- (tumorPurity-1-(bafSites-1)*2^(logRSites/gamma)*((1-tumorPurity)*2+tumorPurity*tumorPloidy))/tumorPurity
  return(cbind(nA,nB))
}

PasteVector <- function(v,sep=""){
  vt <- v[1];
  if(length(v) > 1){
    for(g in 2:length(v)){
      vt <- paste(vt,v[g],sep=sep)
    }
  }
  vt <- paste(vt," EnD",sep="");
  out.v <- sub(" EnD","",vt);
  out.v <- sub("NA , ","",out.v);
  out.v <- sub(" , NA","",out.v);
  out.v <- sub(" , NA , "," , ",out.v);
  return(out.v);
}

getMisMatchPositionsPairwiseAlignment <- function(alignment, chunksize=60, returnlist=FALSE){
  seq1aln <- pattern(alignment) # Get the alignment for the first sequence
  seq2aln <- subject(alignment) # Get the alignment for the second sequence
  #let's check differences across the whole thing
  gapsSeq1        <- countPattern("-",as.character(seq1aln))
  seq1alnresidues <- length(unlist(strsplit(as.character(seq1aln),split="")))-gapsSeq1
  k <- 1
  seq1Positions   <- c()
  for (char in unlist(strsplit(as.character(seq1aln),split="")))
  {
    if(char%in%c('C','G','A','T'))
    {
      seq1Positions <- c(seq1Positions,k)
      k <- k+1
      next;
    }
    if(char%in%c('-'))
    {
      seq1Positions <- c(seq1Positions,k)
      next;
      #
    }
  }
  
  
  k <- 1
  seq2Positions   <- c()
  for (char in unlist(strsplit(as.character(seq2aln),split="")))
  {
    if(char%in%c('C','G','A','T'))
    {
      seq2Positions <- c(seq2Positions,k)
      k <- k+1
      next;
    }
    if(char%in%c('-'))
    {
      seq2Positions <- c(seq2Positions,k)
      next;
      #
    }
  }
  
  diffSeq1 <-  seq1Positions[unlist(strsplit(as.character(seq1aln),split=""))!=unlist(strsplit(as.character(seq2aln),split=""))]
  diffSeq2 <-  seq2Positions[unlist(strsplit(as.character(seq1aln),split=""))!=unlist(strsplit(as.character(seq2aln),split=""))]
  
  diffType1 <- rep(1,length(diffSeq1))
  diffType1[which(unlist(strsplit(as.character(seq1aln),split=""))[unlist(strsplit(as.character(seq1aln),split=""))!=unlist(strsplit(as.character(seq2aln),split=""))]%in%'-')] <- 2
  
  diffType2 <- rep(1,length(diffSeq2))
  diffType2[which(unlist(strsplit(as.character(seq2aln),split=""))[unlist(strsplit(as.character(seq2aln),split=""))!=unlist(strsplit(as.character(seq1aln),split=""))]%in%'-')] <- 2
 
  
  out <- list()
  out$diffSeq1 <- diffSeq1
  out$diffSeq2 <- diffSeq2
  out$diffType1 <- diffType1
  out$diffType2 <- diffType2
  return(out)
}


result <- c()
##start running
HlaFasta   <- read.fasta("test.fasta")
hlaAlleles <-  names(HlaFasta)

for (HLA_gene in c('hla_a','hla_b','hla_c')){
  print(HLA_gene)
  HLA_As <- grep(HLA_gene,hlaAlleles,value=TRUE)
  if (length(HLA_As) < 2){
    next
  }else{
    HLA_type1 <- HLA_As[1]  
    HLA_type2 <- HLA_As[2]
    print(HLA_type1)
    print(HLA_type2)
    HLA_type1Fasta <- HlaFasta[[HLA_type1]]
    HLA_type2Fasta <- HlaFasta[[HLA_type2]]
    # print(HLA_type1Fasta)
    # print(HLA_type2Fasta)
      
    ## 读取.mpileup，得到T和N两个type每个位点的cov
    HLA_type1normal <- read.table(paste(workDir, "/",sampleID,'.',HLA_type1,".n_bwa_sorted_dedup_realign_filt.bam.mpileup",sep=""),sep="\t",stringsAsFactors=FALSE,quote="",fill=TRUE)
    HLA_type1tumor  <- read.table(paste(workDir, "/",sampleID,'.',HLA_type1,".t_bwa_sorted_dedup_realign_filt.bam.mpileup",sep=""),sep="\t",stringsAsFactors=FALSE,quote="",fill=TRUE)
    rownames(HLA_type1normal) <- HLA_type1normal$V2
    rownames(HLA_type1tumor) <- HLA_type1tumor$V2
    HLA_type1normal <- HLA_type1normal[HLA_type1normal$V4>minCoverageFilter,,drop=FALSE] #apply minimum coverage thresholds (we only apply this to the normal for now)
    tmp <- intersect(rownames(HLA_type1tumor),rownames(HLA_type1normal)) 
    HLA_type1tumor  <- HLA_type1tumor[tmp,,drop=FALSE]   
    HLA_type1normalCov <- HLA_type1normal$V4 #rep(0,max(HLA_type1normal$V2,HLA_type1tumor$V2))
    names(HLA_type1normalCov) <- HLA_type1normal$V2 #1:length(HLA_type1tumorCov)
    HLA_type1tumorCov <- rep(0,length(HLA_type1normalCov))      
    names(HLA_type1tumorCov) <- names(HLA_type1normalCov)
    HLA_type1tumorCov[rownames(HLA_type1tumor)] <- HLA_type1tumor$V4
    
    HLA_type2normal <- read.table(paste(workDir, "/",sampleID,'.',HLA_type2,".n_bwa_sorted_dedup_realign_filt.bam.mpileup",sep=""),sep="\t",stringsAsFactors=FALSE,quote="",fill=TRUE)
    HLA_type2tumor  <- read.table(paste(workDir, "/",sampleID,'.',HLA_type2,".t_bwa_sorted_dedup_realign_filt.bam.mpileup",sep=""),sep="\t",stringsAsFactors=FALSE,quote="",fill=TRUE)
    rownames(HLA_type2normal) <- HLA_type2normal$V2
    rownames(HLA_type2tumor) <- HLA_type2tumor$V2
    HLA_type2normal <- HLA_type2normal[HLA_type2normal$V4>minCoverageFilter,,drop=FALSE] #apply minimum coverage thresholds (we only apply this to the normal for now)
    tmp <- intersect(rownames(HLA_type2tumor),rownames(HLA_type2normal))
    HLA_type2tumor  <- HLA_type2tumor[tmp,,drop=FALSE]
    HLA_type2normalCov <- HLA_type2normal$V4 #rep(0,max(HLA_type2normal$V2,HLA_type2tumor$V2))
    names(HLA_type2normalCov) <- HLA_type2normal$V2 #1:length(HLA_type2tumorCov)
    HLA_type2tumorCov <- rep(0,length(HLA_type2normalCov))      
    names(HLA_type2tumorCov) <- names(HLA_type2normalCov)
    HLA_type2tumorCov[rownames(HLA_type2tumor)] <- HLA_type2tumor$V4
    # print("here1")
    # pairmapping得到mismatch position
    seq1 <- PasteVector(toupper(HLA_type1Fasta),sep="")
    seq2 <- PasteVector(toupper(HLA_type2Fasta),sep="")
    sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
    tmp <- pairwiseAlignment(seq1, seq2, substitutionMatrix = sigma, gapOpening = -2,gapExtension = -4, scoreOnly = FALSE,type='local')
    missMatchPositions <- getMisMatchPositionsPairwiseAlignment(tmp,returnlist=TRUE)
    # print("here2")

    # 结果参数
    misMatchCoveredInBoth <- cbind(ifelse(missMatchPositions$diffSeq1%in%names(HLA_type1normalCov),1,0),ifelse(missMatchPositions$diffSeq2%in%names(HLA_type2normalCov),1,0))
    missMatchseq1 <- missMatchPositions$diffSeq1[rowSums(misMatchCoveredInBoth)==2]
    missMatchseq2 <- missMatchPositions$diffSeq2[rowSums(misMatchCoveredInBoth)==2]
    # LossAllele and KeptAllele(所有missmatch位点logR中位数，不使用bin及BAF）：
    HLAtype1Log2MedianCoverageAtSites <- median(log2(c(HLA_type1tumorCov/HLA_type1normalCov)*MultFactor)[names(HLA_type1tumorCov)%in%missMatchseq1], na.rm = TRUE)
    HLAtype2Log2MedianCoverageAtSites <- median(log2(c(HLA_type2tumorCov/HLA_type2normalCov)*MultFactor)[names(HLA_type2tumorCov)%in%missMatchseq2], na.rm = TRUE)
    LossAllele                 <- c(HLA_type1,HLA_type2)[ifelse(HLAtype1Log2MedianCoverageAtSites>HLAtype2Log2MedianCoverageAtSites,2,1)]
    KeptAllele                 <- c(HLA_type1,HLA_type2)[ifelse(HLAtype1Log2MedianCoverageAtSites>HLAtype2Log2MedianCoverageAtSites,1,2)]


    # HLA_type1copyNum_withBAFBin and HLA_type2copyNum_withBAFBin（使用bin和BAF计算missmatch位点logR，取中位数）：
    # binlogRCombined在测试过程由于MultFactor产生微小的误差，误差在alle计算公式中被放大
    startChar <- min(c(as.numeric(names(HLA_type1tumorCov))),as.numeric(names(HLA_type2tumorCov)))
    endChar   <- max(c(as.numeric(names(HLA_type1tumorCov))),as.numeric(names(HLA_type2tumorCov)))
    seqToConsider <- seq(startChar,endChar,by=binSize)
    seqToConsider <- c(seqToConsider[-length(seqToConsider)],endChar+1)   ##表示范围如（1,3000，by=150）
    binLogR       <- c()
    for (i in 1:(length(seqToConsider)-1)){
        PotentialSites   <- as.character(seqToConsider[i]:seqToConsider[i+1])  ##(1,2)
        combinedBinTumor  <- median(as.numeric(as.numeric(HLA_type1tumorCov[names(HLA_type1tumorCov)%in%PotentialSites])), na.rm = TRUE)+median(as.numeric(as.numeric(HLA_type2tumorCov[names(HLA_type2tumorCov)%in%PotentialSites])), na.rm = TRUE)
        combinedBinNormal <- median(as.numeric(as.numeric(HLA_type1normalCov[names(HLA_type1normalCov)%in%PotentialSites])), na.rm= TRUE)+median(as.numeric(as.numeric(HLA_type2normalCov[names(HLA_type2normalCov)%in%PotentialSites])), na.rm = TRUE)
        combinedBinlogR   <- log2(combinedBinTumor/combinedBinNormal*MultFactor) ##测试误差
        type1BinlogR     <- median(log2(as.numeric(as.numeric(HLA_type1tumorCov[names(HLA_type1tumorCov)%in%PotentialSites])/as.numeric(HLA_type1normalCov[names(HLA_type1normalCov)%in%PotentialSites])*MultFactor)), na.rm = TRUE)   
        type2BinlogR     <- median(log2(as.numeric(as.numeric(HLA_type2tumorCov[names(HLA_type2tumorCov)%in%PotentialSites])/as.numeric(HLA_type2normalCov[names(HLA_type2normalCov)%in%PotentialSites])*MultFactor)), na.rm = TRUE)
        binLogR <- rbind(binLogR,cbind(seqToConsider[i],seqToConsider[i+1],combinedBinlogR,type1BinlogR,type2BinlogR))
    }

    tmpOut_cn <- cbind(missMatchseq1
      ,log2(c(HLA_type1tumorCov/HLA_type1normalCov)*MultFactor)[as.character(missMatchseq1)]
      ,HLA_type1tumorCov[as.character(missMatchseq1)]
      ,missMatchseq2
      ,log2(c(HLA_type2tumorCov/HLA_type2normalCov)*MultFactor)[as.character(missMatchseq2)]
      ,HLA_type2tumorCov[as.character(missMatchseq2)]
      ,HLA_type1normalCov[as.character(missMatchseq1)]
      ,HLA_type2normalCov[as.character(missMatchseq2)]
      )
    colnames(tmpOut_cn) <- c('missMatchseq1','logR_type1','TumorCov_type1','missMatchseq2','logR_type2','TumorCov_type2','NormalCov_type1','NormalCov_type2')
    tmpOut_cn  <- tmpOut_cn[!duplicated(tmpOut_cn[,1]),,drop=FALSE]  
    tmpOut_cn  <- tmpOut_cn[!duplicated(tmpOut_cn[,4]),,drop=FALSE]
    combinedTable <- data.frame(tmpOut_cn,stringsAsFactors=FALSE)
    combinedTable$BAFcombined  <- combinedTable$TumorCov_type1/(combinedTable$TumorCov_type1+combinedTable$TumorCov_type2) 

    if(nrow(combinedTable) != 0){
      combinedTable$binlogRCombined <- NA
      combinedTable$binNum          <- NA
      for (i in 1:nrow(combinedTable)){
        combinedTable[i,]$binlogRCombined <- binLogR[which(binLogR[,1]<=as.numeric(combinedTable$missMatchseq1[i])&binLogR[,2]>as.numeric(combinedTable$missMatchseq1[i])),3]
        combinedTable[i,]$binNum       <- binLogR[which(binLogR[,1]<=as.numeric(combinedTable$missMatchseq1[i])&binLogR[,2]>as.numeric(combinedTable$missMatchseq1[i])),1]
      }
    }
    rawValsBin      <- funCalcN_withBAF(combinedTable$binlogRCombined,combinedTable$BAFcombined,tumorPloidy,tumorPurity,gamma)
    combinedTable$nAcombinedBin <- rawValsBin[,1]
    combinedTable$nBcombinedBin <- rawValsBin[,2]
    nA_rawVal_withBAF_bin       <- median(combinedTable[!duplicated(combinedTable$binNum),]$nAcombinedBin, na.rm = TRUE)
    nB_rawVal_withBAF_bin       <- median(combinedTable[!duplicated(combinedTable$binNum),]$nBcombinedBin, na.rm = TRUE)
    HLA_type1copyNum_withBAFBin          <- nA_rawVal_withBAF_bin
    HLA_type2copyNum_withBAFBin          <- nB_rawVal_withBAF_bin



    # Pval_unique (missmatch位点的reads数量计算出的logR的t检验)：
    mismatchPosSeq1 <- HLA_type1tumor[HLA_type1tumor$V2%in%missMatchPositions$diffSeq1,,drop=FALSE]  
    missMatchBed1    <- mismatchPosSeq1[,c(1,2,2)]        ##missMatchBed1是由提取了HLA_type1tumor第一列和两次第二列组成
    missMatchBed1$V2 <- as.numeric(missMatchBed1$V2)-1    ##并对第二列做了处理，变成4列
    missMatchBed1$V3 <- as.numeric(missMatchBed1$V2.1)+1  

    mismatchPosSeq2 <- HLA_type2tumor[HLA_type2tumor$V2%in%missMatchPositions$diffSeq2,,drop=FALSE]
    missMatchBed2    <- mismatchPosSeq2[,c(1,2,2)]
    missMatchBed2$V2 <- as.numeric(missMatchBed2$V2)-1
    missMatchBed2$V3 <- as.numeric(missMatchBed2$V2.1)+1

    write.table(missMatchBed1,file=paste(workDir, '/', sampleID,".",HLA_type1,".bed",sep=""),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
    write.table(missMatchBed2,file=paste(workDir, '/', sampleID,".",HLA_type2,".bed",sep=""),quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)

    ##得到missmatch bed中与bam的交集，写成bed格式
    # use the same mismatch positions, already have the bed to get reads that do overlap a mismatch
    Type1NormalCmd <- paste("bedtools intersect -loj -bed -b ",workDir, '/', sampleID,'.',HLA_type1,".n_bwa_sorted_dedup_realign_filt.bam"," -a ",workDir, '/', sampleID,".",HLA_type1,".bed"," > ",workDir, '/', sampleID,".",HLA_type1,".normal.mismatch.reads.bed",sep="")
    Type1TumorCmd <- paste("bedtools intersect -loj -bed -b ",workDir, '/', sampleID,'.',HLA_type1,".t_bwa_sorted_dedup_realign_filt.bam"," -a ",workDir, '/', sampleID,".",HLA_type1,".bed"," > ",workDir, '/', sampleID,".",HLA_type1,".tumor.mismatch.reads.bed",sep="")
    Type2NormalCmd <- paste("bedtools intersect -loj -bed -b ",workDir, '/', sampleID,'.',HLA_type2,".n_bwa_sorted_dedup_realign_filt.bam"," -a ",workDir, '/', sampleID,".",HLA_type2,".bed"," > ",workDir, '/', sampleID,".",HLA_type2,".normal.mismatch.reads.bed",sep="")
    Type2TumorCmd <- paste("bedtools intersect -loj -bed -b ",workDir, '/', sampleID,'.',HLA_type2,".t_bwa_sorted_dedup_realign_filt.bam"," -a ",workDir, '/', sampleID,".",HLA_type2,".bed"," > ",workDir, '/', sampleID,".",HLA_type2,".tumor.mismatch.reads.bed",sep="")
    system(Type1NormalCmd)
    system(Type1TumorCmd)
    system(Type2NormalCmd)
    system(Type2TumorCmd)

  # only take unique reads
    for (i in list.files(pattern = 'mismatch.reads.bed', full.names = TRUE)){
      if(file.size(i) > 0){
      x <- read.csv(i, sep = '\t', as.is = TRUE, header = FALSE)
      x <- x[!is.na(x$V2),]
      x <- x[!duplicated(x$V8),]
      } else{
        x <- t(c(paste('V', seq(1, 10), sep = '')))
      }
      write.table(x, file = gsub(pattern = 'mismatch.reads.bed', replacement = 'mismatch.unique.reads.bed', x = i), sep = '\t', quote = FALSE, row.names = FALSE)
      # system(paste("mv ",i," ",workDir,"./MismatchReadsBbed",sep="")) 
    }
    
    HLA_type1normal_unique             <- read.table(paste(workDir, '/', sampleID, '.', HLA_type1, '.normal.mismatch.unique.reads.bed', sep = ''), sep= '\t', header = TRUE, as.is = TRUE, comment.char = '')
    HLA_type1normalCov_mismatch        <- HLA_type1normalCov[names(HLA_type1normalCov)%in%missMatchPositions$diffSeq1]
    HLA_type1normalCov_mismatch_unique <- HLA_type1normalCov_mismatch
    HLA_type1normalCov_mismatch_unique <- sapply(X = names(HLA_type1normalCov_mismatch_unique), FUN = function(x) {return(table(HLA_type1normal_unique$V3)[x])} )  
    names(HLA_type1normalCov_mismatch_unique) <- names(HLA_type1normalCov_mismatch)
    HLA_type1normalCov_mismatch_unique[is.na(HLA_type1normalCov_mismatch_unique)] <- 0
    if(length(HLA_type1normalCov_mismatch_unique) != 0){
      HLA_type1normalCov_mismatch_unique <- HLA_type1normalCov_mismatch_unique + 1
    }

    HLA_type1tumor_unique             <- read.table(paste(workDir, '/', sampleID, '.', HLA_type1, '.tumor.mismatch.unique.reads.bed', sep = ''), sep= '\t', header = TRUE, as.is = TRUE, comment.char = '')
    HLA_type1tumorCov_mismatch        <- HLA_type1tumorCov[names(HLA_type1tumorCov)%in%missMatchPositions$diffSeq1]
    HLA_type1tumorCov_mismatch_unique <- HLA_type1tumorCov_mismatch
    HLA_type1tumorCov_mismatch_unique <- sapply(X = names(HLA_type1tumorCov_mismatch_unique), FUN = function(x) {return(table(HLA_type1tumor_unique$V3)[x])} )  
    names(HLA_type1tumorCov_mismatch_unique) <- names(HLA_type1tumorCov_mismatch)
    HLA_type1tumorCov_mismatch_unique[is.na(HLA_type1tumorCov_mismatch_unique)] <- 0
    if(length(HLA_type1tumorCov_mismatch_unique) != 0){
      HLA_type1tumorCov_mismatch_unique <- HLA_type1tumorCov_mismatch_unique + 1
    }

    HLA_type2normal_unique             <- read.table(paste(workDir, '/', sampleID, '.', HLA_type2, '.normal.mismatch.unique.reads.bed', sep = ''), sep= '\t', header = TRUE, as.is = TRUE, comment.char = '')
    HLA_type2normalCov_mismatch        <- HLA_type2normalCov[names(HLA_type2normalCov)%in%missMatchPositions$diffSeq2]
    HLA_type2normalCov_mismatch_unique <- HLA_type2normalCov_mismatch
    HLA_type2normalCov_mismatch_unique <- sapply(X = names(HLA_type2normalCov_mismatch_unique), FUN = function(x) {return(table(HLA_type2normal_unique$V3)[x])} )  
    names(HLA_type2normalCov_mismatch_unique) <- names(HLA_type2normalCov_mismatch)
    HLA_type2normalCov_mismatch_unique[is.na(HLA_type2normalCov_mismatch_unique)] <- 0
    if(length(HLA_type2normalCov_mismatch_unique) != 0){
      HLA_type2normalCov_mismatch_unique <- HLA_type2normalCov_mismatch_unique + 1
    }

    HLA_type2tumor_unique             <- read.table(paste(workDir, '/', sampleID, '.', HLA_type2, '.tumor.mismatch.unique.reads.bed', sep = ''), sep= '\t', header = TRUE, as.is = TRUE, comment.char = '')
    HLA_type2tumorCov_mismatch        <- HLA_type2tumorCov[names(HLA_type2tumorCov)%in%missMatchPositions$diffSeq2]
    HLA_type2tumorCov_mismatch_unique <- HLA_type2tumorCov_mismatch
    HLA_type2tumorCov_mismatch_unique <- sapply(X = names(HLA_type2tumorCov_mismatch_unique), FUN = function(x) {return(table(HLA_type2tumor_unique$V3)[x])} )  
    names(HLA_type2tumorCov_mismatch_unique) <- names(HLA_type2tumorCov_mismatch)
    HLA_type2tumorCov_mismatch_unique[is.na(HLA_type2tumorCov_mismatch_unique)] <- 0
    if(length(HLA_type2tumorCov_mismatch_unique) != 0){
      HLA_type2tumorCov_mismatch_unique <- HLA_type2tumorCov_mismatch_unique + 1
    }

    tmpOut_unique <- cbind(missMatchseq1,log2(c(HLA_type1tumorCov_mismatch_unique/HLA_type1normalCov_mismatch_unique + 0.0001)*MultFactor)[as.character(missMatchseq1)],missMatchseq2,log2(c(HLA_type2tumorCov_mismatch_unique/HLA_type2normalCov_mismatch_unique+0.0001)*MultFactor)[as.character(missMatchseq2)])
    tmpOut_unique  <- tmpOut_unique[!duplicated(tmpOut_unique[,1]),,drop=FALSE]
    tmpOut_unique  <- tmpOut_unique[!duplicated(tmpOut_unique[,3]),,drop=FALSE]
    if(nrow(tmpOut_unique) > 1){
      PairedTtest_unique <- t.test(tmpOut_unique[,2],tmpOut_unique[,4],paired=TRUE)
    }
    PVal_unique  <- PairedTtest_unique$p.value

    ##单个HLA结果表格
    HLAout <- paste(workDir, '/', sampleID,'.',HLA_gene,'.',minCoverageFilter,".DNA.IntegerCPN_CI.xls",sep="")
    write.table(combinedTable,file=HLAout,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

    PatientOutPut <- cbind(HLA_type1,HLA_type2,PVal_unique,HLAtype1Log2MedianCoverageAtSites,HLAtype2Log2MedianCoverageAtSites,LossAllele,KeptAllele,HLA_type1copyNum_withBAFBin,HLA_type2copyNum_withBAFBin)
    result <- rbind(result,PatientOutPut)
  }
}

##结果汇总
HLAoutall <- paste(workDir, '/', sampleID,'.',minCoverageFilter,".DNA.HLAlossPrediction_CI.xls",sep="")
write.table(result,file=HLAoutall,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)      
