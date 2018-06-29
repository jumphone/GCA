library(gplots)
a=as.matrix(read.table('MODE_MAT.txt',header=T,row.names=1))
stem_score=read.table('stem_score.txt')[,2]
ac_score=read.table('AC_score.txt')[,2]
oc_score=read.table('OC_score.txt')[,2]


SUM=apply(a,2,sum)
B= which(SUM>=3)
b=a[,B]

stem_score_b=stem_score[B]
ac_score_b=ac_score[B]
oc_score_b=oc_score[B]

length(stem_score_b)

RSUM=apply(b,1,sum)
b=b[which(RSUM>0),]
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

PCAPlot(object = EXP, dim.1 = 1, dim.2 = 2)

PCUSE=1:10
EXP = RunTSNE(object = EXP, dims.use = PCUSE, do.fast = TRUE,check_duplicates = FALSE )

RES=0.6
EXP <- FindClusters(object = EXP, reduction.type = "pca", dims.use = PCUSE,  resolution = RES, print.output = 0, save.SNN = TRUE,force.recalc =T)

TSNEPlot(object = EXP,do.label=T)

#VlnPlot(object = EXP, features.plot = c('stem.score'),do.sort=T)



pdf('STEM_AC_OC.pdf',width=10,height=10)
COMBINE=c()
##########STEM###########
stem_score_b=stem_score[B]
med_stem_score=c()
wilcox_test_p=c()
i=0
while(i<length(table(EXP@ident))){
med_stem_score=c(med_stem_score, median(stem_score_b[which(EXP@ident==i)]))
wilcox_test_p=c(wilcox_test_p,wilcox.test(stem_score_b[which(EXP@ident==i)],stem_score_b)$p.value ) 
i=i+1
}
neg_log_2_p=-log(wilcox_test_p,2)
neg_log_2_adjp= -log(p.adjust(wilcox_test_p,method='fdr'),2)
direction= (med_stem_score - median(stem_score_b))/abs((med_stem_score - median(stem_score_b)))
#plot(neg_log_2_p*direction,pch=16)
#abline(h=0)
#abline(h=-log(0.05,2),lty=3)
#abline(h=log(0.05,2),lty=3)
COL=rep('black',length(neg_log_2_adjp))
COL[which(neg_log_2_adjp > -log(0.05,2))]='red'
plot(x=c(0:12),y=neg_log_2_adjp*direction,pch=16,ylim=c(-150,50),col=COL,cex=3,main='STEM')
abline(h=0)
abline(h=-log(0.05,2),lty=3)
abline(h=log(0.05,2),lty=3)

COMBINE=cbind(COMBINE,neg_log_2_adjp*direction)
#################################

##########AC###########
stem_score_b=ac_score_b
med_stem_score=c()
wilcox_test_p=c()
i=0
while(i<length(table(EXP@ident))){
med_stem_score=c(med_stem_score, median(stem_score_b[which(EXP@ident==i)]))
wilcox_test_p=c(wilcox_test_p,wilcox.test(stem_score_b[which(EXP@ident==i)],stem_score_b)$p.value ) 
i=i+1
}
neg_log_2_p=-log(wilcox_test_p,2)
neg_log_2_adjp= -log(p.adjust(wilcox_test_p,method='fdr'),2)
direction= (med_stem_score - median(stem_score_b))/abs((med_stem_score - median(stem_score_b)))
#plot(neg_log_2_p*direction,pch=16)
#abline(h=0)
#abline(h=-log(0.05,2),lty=3)
#abline(h=log(0.05,2),lty=3)
COL=rep('black',length(neg_log_2_adjp))
COL[which(neg_log_2_adjp > -log(0.05,2))]='red'
plot(x=c(0:12),y=neg_log_2_adjp*direction,pch=16,col=COL,cex=3,main='AC')
abline(h=0)
abline(h=-log(0.05,2),lty=3)
abline(h=log(0.05,2),lty=3)
COMBINE=cbind(COMBINE,neg_log_2_adjp*direction)
#################################

##########OC###########
stem_score_b=oc_score_b
med_stem_score=c()
wilcox_test_p=c()
i=0
while(i<length(table(EXP@ident))){
med_stem_score=c(med_stem_score, median(stem_score_b[which(EXP@ident==i)]))
wilcox_test_p=c(wilcox_test_p,wilcox.test(stem_score_b[which(EXP@ident==i)],stem_score_b)$p.value ) 
i=i+1
}
neg_log_2_p=-log(wilcox_test_p,2)
neg_log_2_adjp= -log(p.adjust(wilcox_test_p,method='fdr'),2)
direction= (med_stem_score - median(stem_score_b))/abs((med_stem_score - median(stem_score_b)))
#plot(neg_log_2_p*direction,pch=16)
#abline(h=0)
#abline(h=-log(0.05,2),lty=3)
#abline(h=log(0.05,2),lty=3)
COL=rep('black',length(neg_log_2_adjp))
COL[which(neg_log_2_adjp > -log(0.05,2))]='red'
plot(x=c(0:12),y=neg_log_2_adjp*direction,pch=16,col=COL,cex=3,main='OC')
abline(h=0)
abline(h=-log(0.05,2),lty=3)
abline(h=log(0.05,2),lty=3)
COMBINE=cbind(COMBINE,neg_log_2_adjp*direction)
#################################

NEW_COMBINE=COMBINE
NEW_COMBINE[which(COMBINE> -log(0.05,2))]=1
NEW_COMBINE[which(COMBINE > log(0.05,2) & COMBINE < -log(0.05,2))]=0
NEW_COMBINE[which(COMBINE< log(0.05,2))]= 0
#heatmap.2(NEW_COMBINE,scale='none')


rownames(NEW_COMBINE)=as.character(c(0:12))
colnames(NEW_COMBINE)=c('STEM','AC','OC')
library('gplots')
#heatmap.2(NEW_COMBINE,scale='none',trace='none',col=colorRampPalette(c('blue','grey80','red')),cexCol=1)
NEW_NEW_COMBINE=NEW_COMBINE[,c(2,1,3)]
NEW_NEW_COMBINE=NEW_NEW_COMBINE[c(5,7,11,4,1,3,0,2,6,8,10,9,12)+1,]

#colnames(NEW_NEW_COMBINE)=c('AC','STEM','OC')
heatmap.2(NEW_NEW_COMBINE,scale='none',trace='none',col=colorRampPalette(c('grey80','red')),cexCol=1,Colv=F,Rowv=F,dendrogram='none')

dev.off()


#################################

TSNEPlot(object = EXP,do.label=T,label.size=10)




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


plot(neg_log_2_adjp*direction,pch=16,ylim=c(-150,50),col=COL)
abline(h=0)
abline(h=-log(0.05,2),lty=3)
abline(h=log(0.05,2),lty=3)

dev.off()

write.table(top10,file='top10.txt',row.names=T,col.names=T,quote=F,sep='\t')


