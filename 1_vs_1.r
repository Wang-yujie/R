# in R command prompt
##Rscript $0 Control1.txt Control2.txt Control3.txt Treat1.txt treat2.txt treat3.txt
args<-commandArgs(T)
Control=read.table(args[1],sep = '\t', row.names=1)
Treat=read.table(args[2],sep = '\t', row.names=1)
Treat_vs_Control=cbind(Treat, Control)
colnames(Treat_vs_Control)=c("Treat", "Control")
write.table(Treat_vs_Control,file="Treat_vs_Control.test", quote = F,sep = '\t')
countdata=read.table("Treat_vs_Control.test",header=T,row.names=1)
condition = c(rep('Treat',1), rep('Control',1))
library(DESeq2)
coldata <- data.frame(row.names = colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~condition)
dds <- DESeq(dds)
res <- results(dds,contrast = c("condition", "Treat", "Control"))
res <- res[order(res$pvalue),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
rownames(resdata)<-resdata[,1]
resdata<-resdata[,2:dim(resdata)[2]]
diff<-resdata[! is.na(resdata[,5]),]
diff<-diff[diff[,5]<0.05,]
up<-diff[diff[,4]>0,]
down<-diff[diff[,4]<0,]
write.table(rownames(down),file="gene_down.test",quote=F,row.names = F,col.names = F)
write.table(rownames(up),file="gene_up.test",quote=F,row.names = F,col.names = F)
newresdata=data.frame(gene_id=rownames(resdata),resdata)
write.table(newresdata, file="Treat_vs_Control_diff.test",sep="\t",quote=F,row.names = FALSE)