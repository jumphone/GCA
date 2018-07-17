library(Seurat)
library(dplyr)
library(Matrix)

exp_data=Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
EXP = CreateSeuratObject(raw.data = exp_data, min.cells = 3, min.genes=100)
dim(EXP@data)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = EXP@data), value = TRUE)
percent.mito <- colSums(EXP@data[mito.genes, ]) / colSums(EXP@data)
EXP <- AddMetaData(object = EXP, metadata = percent.mito, col.name = "percent.mito")

pdf('Seurat_QC.pdf',width=30,height=15)
VlnPlot(object = EXP, features.plot = c("nGene", "nUMI",'percent.mito'), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = EXP, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = EXP, gene1 = "nUMI", gene2 = "nGene")
dev.off()

EXP=FilterCells(object = EXP, subset.names = c("nGene", "percent.mito"), low.thresholds = c(100, -Inf), high.thresholds = c(5000, 0.05))
dim(EXP@data)
EXP <- NormalizeData(object =EXP, normalization.method = "LogNormalize",  scale.factor = 10000)
all_gene=rownames(EXP@data)
EXP <- ScaleData(object = EXP, vars.to.regress = c("nUMI", "percent.mito"), genes.use =all_gene)
#######DeScale##############
MIN=apply(EXP@scale.data,1,min)
OUT=as.matrix(EXP@scale.data)
i=1
while(i <= length(all_gene)){
    OUT[i,]=as.matrix(EXP@scale.data)[i,]-MIN[i]
    i=i+1
}
OUT[which(EXP@data==0)]=0
#############
summary(OUT[,1])
write.table(file='normalized_exp.txt',OUT,row.names=T,col.names=T,quote=F,sep='\t')
save(EXP, file = "Seurat_EXP.Robj")
