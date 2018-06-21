library(gplots)
a=as.matrix(read.table('MODE_MAT.txt',header=T,row.names=1))
stem_score=read.table('')

SUM=apply(a,2,sum)
B= which(SUM>=3)
b=a[,B]
RSUM=apply(b,1,sum)
b=b[which(RSUM>0),]
heatmap(b,Rowv=T,Colv=F,scale='row',labCol='',margins=c(10,10))
all_gene=rownames(b)

library(Seurat)

EXP = CreateSeuratObject(raw.data = b, min.cells = 0, min.genes=0)

EXP=NormalizeData(object = EXP, normalization.method = "LogNormalize", scale.factor = 10000)

EXP = ScaleData(object = EXP,, genes.use = all_gene)

PCNUM=20
EXP <- RunPCA(object = EXP, pc.genes = all_gene, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )
PCElbowPlot(object = EXP,num.pc=PCNUM)

PCUSE=1:3
EXP = RunTSNE(object = EXP, dims.use = PCUSE, do.fast = TRUE,check_duplicates = FALSE )

RES=0.3
EXP <- FindClusters(object = EXP, reduction.type = "pca", dims.use = PCUSE,  resolution = RES, print.output = 0, save.SNN = TRUE,force.recalc =T)

TSNEPlot(object = EXP,do.label=T)


pbmc=EXP
library(dplyr)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.1, thresh.use = 0.2)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,col.low = "grey90", col.mid = "grey75", col.high = "red" )
