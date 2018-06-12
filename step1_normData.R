library(Seurat)
library(dplyr)
library(Matrix)

#args = commandArgs(trailingOnly=TRUE)

exp_data=read.table('10X_Tumor72017_matrix.txt.uniq.txt',header=T,row.names=1,sep='\t')


EXP = CreateSeuratObject(raw.data = exp_data, min.cells = 3, min.genes=200)
mito.genes <- grep(pattern = "^mt-", x = rownames(x = EXP@data), value = TRUE)
percent.mito <- colSums(EXP@data[mito.genes, ]) / colSums(EXP@data)
EXP <- AddMetaData(object = EXP, metadata = percent.mito, col.name = "percent.mito")

pdf('Seurat_QC.pdf',width=30,height=15)
VlnPlot(object = EXP, features.plot = c("nGene", "nUMI",'percent.mito'), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = EXP, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = EXP, gene1 = "nUMI", gene2 = "nGene")
dev.off()

EXP=FilterCells(object = EXP, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(6000, 0.05))
EXP=NormalizeData(object = EXP, normalization.method = "LogNormalize", scale.factor = 10000)

save(EXP, file = "Seurat_EXP.Robj")
OUT=as.matrix(EXP@data)
write.table(file='normalized_exp.txt',OUT,row.names=T,col.names=T,quote=F,sep='\t')
