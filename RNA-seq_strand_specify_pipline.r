# in R command prompt
##Rscript $0 Control1.txt Control2.txt Control3.txt Treat1.txt treat2.txt treat3.txt
args<-commandArgs(T)
Control1=read.table(args[1],sep = '\t', row.names=1)
Control2=read.table(args[2],sep = '\t', row.names=1)
Control3=read.table(args[3],sep = '\t', row.names=1)
Control=cbind(Control1, Control2, Control3)
colnames(Control)=c("Control1", "Control2", "Control3")
Treat1=read.table(args[4],sep = '\t', row.names=1)
Treat2=read.table(args[5],sep = '\t', row.names=1)
Treat3=read.table(args[6],sep = '\t', row.names=1)
Treat=cbind(Treat1, Treat2, Treat3)
colnames(Treat)=c("Treat1", "Treat2" ,"Treat3")
Treat_vs_Control=cbind(Treat, Control)
write.table(Treat_vs_Control,file="Treat_vs_Control.txt", quote = F,sep = '\t')
countdata=read.table("Treat_vs_Control.txt",header=T,row.names=1)
condition = c(rep('Treat',3), rep('Control',3))
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
write.table(rownames(down),file="gene_down.txt",quote=F,row.names = F,col.names = F)
write.table(rownames(up),file="gene_up.txt",quote=F,row.names = F,col.names = F)
newresdata=data.frame(gene_id=rownames(resdata),resdata)
write.table(newresdata, file="Treat_vs_Control_diff.txt",sep="\t",quote=F,row.names = FALSE)