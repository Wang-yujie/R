args<-commandArgs(T)
library(pheatmap)
pdata=read.table(args[1],header=T,row.names=1)
pdata<-pdata[pdata$pvalue<0.05,]
pdata<-pdata[! is.na(pdata[,5]),]
pdata$ave <- apply(pdata[,7:12],1,mean)
pdata[,7:12] <- pdata[,7:12]/pdata$ave
km <- kmeans(pdata[,7:12],2)
pdata$cl <- km$cluster
pdata <- pdata[order(pdata$cl),]
ann <- data.frame(cluster=pdata$cl)
rownames(ann) <- rownames(pdata)
ann_col=data.frame(sample=c(rep('control',3),rep('treat',3)))
rownames(ann_col)=colnames(pdata)[7:12]
pdf("args[1].pdf")
pheatmap(pdata[,7:12],col=colorRampPalette(c('blue','white','red'))(50),scale='none',cluster_rows=F,cluster_cols=F,show_rownames=F,cellwidth=30,annotation_row=ann,annotation_col=ann_col)
dev.off()