args<-commandArgs(T)
library(pathview)
library(Category)
library(org.Hs.eg.db)
library(GOstats)
library(KEGG.db)
library(GO.db)

data=read.table(args[1], sep="\t")
a=t(data)
b=as.vector(a)
genes=b
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
goAnn <- get("org.Hs.egGO")
universe <- Lkeys(goAnn)
params <- new("GOHyperGParams",
              geneIds=entrezIDs,
              universeGeneIds=universe,
              annotation="org.Hs.eg.db",
              ontology="BP",
              pvalueCutoff=0.05,
              testDirection="over")
over <- hyperGTest(params)
bp <- summary(over)
glist <- geneIdsByCategory(over)
glist <- sapply(glist, function(.ids) {
  .sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
  .sym[is.na(.sym)] <- .ids[is.na(.sym)]
  paste(.sym, collapse=";")
})
bp$Symbols <- glist[as.character(bp$GOBPID)]
write.table(bp,file="go.xls",col.names=T,row.names=F,sep="\t")  #写出GO表格，画图要用

##KEGG
keggAnn <- get("org.Hs.egPATH")
universe <- Lkeys(keggAnn)
params <- new("KEGGHyperGParams", geneIds=entrezIDs, universeGeneIds=universe, annotation="org.Hs.eg.db", categoryName="KEGG", pvalueCutoff=0.05,testDirection="over")
over <- hyperGTest(params)
kegg <- summary(over)
glist <- geneIdsByCategory(over)
glist <- sapply(glist, function(.ids) {
  .sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA)
 	.sym[is.na(.sym)] <- .ids[is.na(.sym)]
 	paste(.sym, collapse=";")
})
kegg$Symbols <- glist[as.character(kegg$KEGGID)]
write.table(kegg,file="kegg.xls",col.names=T,row.names=F,sep="\t") #写出KEGG表格

gIds <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
gEns <- unlist(gIds)
gene.data <- rep(1, length(gEns))
names(gene.data) <- gEns
#自动从KEGG上下载各种pathway图和表，一般没人看这玩意儿
for(i in 1:nrow(kegg)){
  pv.out <- pathview(gene.data, pathway.id=as.character(kegg$KEGGID)[i], species="hsa",
                     kegg.native=T,low = list(gene = "green"),
                     mid =list(gene = "gray"), high = list(gene = "green"),plot.col.key=FALSE)
}

#使用ggplot2软件画GO和KEGG直方图了
library(ggplot2)
bp = read.table("go.xls", header=T,row.names=1, sep="\t")
bp_sort=bp[order(bp$Pvalue),]
tmp=bp_sort[1:20,] ## plot top 20 terms
tmp=tmp[order(tmp$Pvalue,decreasing = T),]
tmp$Term=factor(tmp$Term,levels = tmp$Term)
ggplot(tmp)+geom_bar(aes(x=Term, y=-log(Pvalue,10)), stat="identity")+theme(axis.text.x=element_text(angle = 90,hjust=1))+coord_flip()
goplotname="GOplot.pdf"
ggsave(width = 8,height = 8,filename =goplotname)

kg = read.table("kegg.xls", header=T,row.names=1, sep="\t")
kg_sort=kg[order(kg$Pvalue),]
tmp=kg_sort[1:20,] ## plot top 20 terms
tmp=tmp[order(tmp$Pvalue,decreasing = T),]
tmp$Term=factor(tmp$Term,levels = tmp$Term)
ggplot(tmp)+geom_bar(aes(x=Term, y=-log(Pvalue,10)), stat="identity")+theme(axis.text.x=element_text(angle = 90,hjust=1))+coord_flip()
goplotname="KEGGplot.pdf"
ggsave(width = 10,height = 10,filename =goplotname)