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
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,col.low = "grey90", col.mid = "grey75", col.high = "red",cex.row=6 )





