library(gplots)
a=as.matrix(read.table('MODE_MAT.txt',header=T,row.names=1))
stem_score=read.table('stem_score.txt')[,2]

SUM=apply(a,2,sum)
B= which(SUM>=3)
b=a[,B]
stem_score_b=stem_score[B]
length(stem_score_b)

RSUM=apply(b,1,sum)
b=b[which(RSUM>0),]
#heatmap(b,Rowv=T,Colv=F,scale='row',labCol='',margins=c(10,10))
all_gene=rownames(b)

library(Seurat)

EXP = CreateSeuratObject(raw.data = b, min.cells = 0, min.genes=0)

EXP <- AddMetaData(object = EXP, metadata = stem_score_b, col.name = "stem.score")
EXP@meta.data$stem.score=stem_score_b

EXP=NormalizeData(object = EXP, normalization.method = "LogNormalize", scale.factor = 10000)

EXP = ScaleData(object = EXP,, genes.use = all_gene)

PCNUM=40
EXP <- RunPCA(object = EXP, pc.genes = all_gene, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )
PCElbowPlot(object = EXP,num.pc=PCNUM)

#EXP <- JackStraw(object = EXP, num.replicate = 100, display.progress = FALSE)
#JackStrawPlot(object = EXP, PCs = 1:PCNUM)

PCAPlot(object = EXP, dim.1 = 1, dim.2 = 2)

PCUSE=1:10
EXP = RunTSNE(object = EXP, dims.use = PCUSE, do.fast = TRUE,check_duplicates = FALSE )

RES=0.6
EXP <- FindClusters(object = EXP, reduction.type = "pca", dims.use = PCUSE,  resolution = RES, print.output = 0, save.SNN = TRUE,force.recalc =T)

TSNEPlot(object = EXP,do.label=T)

VlnPlot(object = EXP, features.plot = c('stem.score'))


med_stem_score=c()
i=0
while(i<length(table(EXP@ident))){

med_stem_score=c(med_stem_score, median(stem_score_b[which(EXP@ident==i)]))
i=i+1
}

O=order(med_stem_score,decreasing=T)
plot(med_stem_score[O])

boxplot( 
  stem_score_b[which(EXP@ident==4-1)], 
  stem_score_b[which(EXP@ident==5-1)], 
  stem_score_b[which(EXP@ident==2-1)], 
  stem_score_b[which(EXP@ident==13-1)], 
  stem_score_b[which(EXP@ident==8-1)], 
  stem_score_b[which(EXP@ident==3-1)], 
  stem_score_b[which(EXP@ident==6-1)], 
  stem_score_b[which(EXP@ident==9-1)], 
  stem_score_b[which(EXP@ident==10-1)], 
  stem_score_b[which(EXP@ident==1-1)], 
  stem_score_b[which(EXP@ident==7-1)], 
  stem_score_b[which(EXP@ident==11-1)], 
  stem_score_b[which(EXP@ident==12-1)], 
    outline=F   )




EXP@scale.data=as.matrix(EXP@data)
pbmc=EXP
library(dplyr)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.05, thresh.use = 0.1)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,col.low = "grey90", col.mid = "grey75", col.high = "red",cex.row=6 )







