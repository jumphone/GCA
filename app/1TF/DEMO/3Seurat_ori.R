

################################################################################
################################################################################
################################################################################
PCNUM=40
PCUSE=1:10
RES=0.6

exp_data=as.matrix(read.table('../GSE70630_OG_processed_data_v2.txt.cleaned.txt',header=T,row.names=1))
tf_ident=read.table('IDENT.txt')
new_exp_data=exp_data[,which(colnames(exp_data) %in% tf_ident[,1])]
new_exp_data[is.na(new_exp_data)]=0
EXP = CreateSeuratObject(raw.data = new_exp_data, min.cells = -1, min.genes=-1)
EXP=NormalizeData(object = EXP, normalization.method = "LogNormalize", scale.factor = 10000)
pdf('EXP_VarGene.pdf')
EXP <- FindVariableGenes(object = EXP, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff =0, y.cutoff = 0.8)
dev.off()
length(x=EXP@var.genes)

EXP = ScaleData(object = EXP, genes.use = EXP@var.genes)
EXP <- RunPCA(object = EXP, pc.genes = EXP@var.genes, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )

EXP <- FindClusters(object = EXP, reduction.type = "pca", dims.use = PCUSE,  resolution = RES, print.output = 0, save.SNN = TRUE)
EXP <- RunTSNE(object = EXP, dims.use = PCUSE, do.fast = TRUE)
pdf('EXP_TSNE.pdf')
TSNEPlot(object = EXP,do.label=T)
dev.off()
write.table(file='EXP_IDENT.txt', EXP@ident,row.names=T,col.names=F,sep='\t',quote=F)

TAG=read.table('MODE_MAT.txt.tag')[,1]
new_EXP=EXP
TAG=as.factor(TAG)
names(TAG)=names(new_EXP@ident)
new_EXP@ident=TAG
pdf('EXP_BATCH.pdf',width=7,height=7)
TSNEPlot(object = new_EXP,do.label=F)
dev.off()




################################################################################
################################################################################
################################################################################


library(gplots)
library(Seurat)

a=as.matrix(read.table('MODE_MAT.txt',header=T,row.names=1))
RSUM=apply(b,1,sum)
b=b[which(RSUM>0),]
all_gene=rownames(b)

exp_ident=read.table('EXP_IDENT.txt')

new_exp_data=b[,which(colnames(b) %in% exp_ident[,1])]

EXP = CreateSeuratObject(raw.data = new_exp_data, min.cells = 0, min.genes=0)


tmp_ident=exp_ident[,2]
tmp_ident=as.factor(tmp_ident)
names(tmp_ident)=names(EXP@ident)
EXP@ident=tmp_ident
head(EXP@ident)
EXP = ScaleData(object = EXP, genes.use = all_gene)


pbmc=EXP
library(dplyr)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0, thresh.use = 1)

pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

write.table(file='EXP_MARKER.txt', pbmc.markers,row.names=T,col.names=T,quote=F,sep='\t')

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

all_gene=rownames(EXP@data)



pdf('EXPC_TF_MARKER.pdf',width=7,height=7)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,col.low = "grey90", col.mid = "grey75", col.high = "red",cex.row=6,cex.col=1)
dev.off()




###########
###########
###########
###########
###########
library(gplots)
library(Seurat)

exp_data=as.matrix(read.table('GSE70630_OG_processed_data_v2.txt.cleaned.txt',header=T,row.names=1))
tf_ident=read.table('IDENT.txt')

new_exp_data=exp_data[,which(colnames(exp_data) %in% tf_ident[,1])]

EXP = CreateSeuratObject(raw.data = new_exp_data, min.cells = 0, min.genes=0)


tmp_ident=tf_ident[,2]
tmp_ident=as.factor(tmp_ident)
names(tmp_ident)=names(EXP@ident)
EXP@ident=tmp_ident
head(EXP@ident)



pbmc=EXP
library(dplyr)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.1)

pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

write.table(file='EXP_MARKER.txt', pbmc.markers,row.names=T,col.names=T,quote=F,sep='\t')

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

all_gene=rownames(EXP@data)
EXP = ScaleData(object = EXP, genes.use = all_gene)


pdf('GENE_MARKER.pdf',width=15,height=15)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,col.low = "grey90", col.mid = "grey75", col.high = "red",cex.row=6 )
dev.off()







