library(gplots)
a=as.matrix(read.table('MODE_MAT.txt',header=T,row.names=1))


SUM=apply(a,2,sum)
B= which(SUM>=3)
b=a[,B]



RSUM=apply(b,1,sum)
b=b[which(RSUM>0),]
all_gene=rownames(b)

library(Seurat)

EXP = CreateSeuratObject(raw.data = b, min.cells = 0, min.genes=0)

EXP=NormalizeData(object = EXP, normalization.method = "LogNormalize", scale.factor = 10000)

EXP = ScaleData(object = EXP,, genes.use = all_gene)

PCNUM=40
EXP <- RunPCA(object = EXP, pc.genes = all_gene, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )
PCElbowPlot(object = EXP,num.pc=PCNUM)

PCAPlot(object = EXP, dim.1 = 1, dim.2 = 2)

PCUSE=1:10
EXP = RunTSNE(object = EXP, dims.use = PCUSE, do.fast = TRUE,check_duplicates = FALSE )

RES=0.6
EXP <- FindClusters(object = EXP, reduction.type = "pca", dims.use = PCUSE,  resolution = RES, print.output = 0, save.SNN = TRUE,force.recalc =T)

TSNEPlot(object = EXP,do.label=T)




EXP@scale.data=as.matrix(EXP@data)
pbmc=EXP
library(dplyr)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.05, thresh.use = 0.1)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,col.low = "grey90", col.mid = "grey75", col.high = "red",cex.row=6 )



pdf('OUTPUT.pdf',width=15,height=15)

TSNEPlot(object = EXP,do.label=T)
VlnPlot(object = EXP, features.plot = c('stem.score'),do.sort=T)
plot(med_stem_score[O],pch=16)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,col.low = "grey90", col.mid = "grey90", col.high = "red",cex.row=6 )
dev.off()

write.table(top10,file='top10.txt',row.names=T,col.names=T,quote=F,sep='\t')


