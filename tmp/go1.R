library(SingleCellExperiment)
library(SC3)
library(scater)


a=as.matrix(read.table('MODE_MAT.txt',header=T,row.names=1))
stem_score=read.table('stem_score.txt')[,2]

SUM=apply(a,2,sum)
B= which(SUM>=2)
b=a[,B]
stem_score_b=stem_score[B]
length(stem_score_b)

RSUM=apply(b,1,sum)
b=b[which(RSUM>0),]
#heatmap(b,Rowv=T,Colv=F,scale='row',labCol='',margins=c(10,10))
all_gene=rownames(b)

COLDATA=as.matrix(stem_score_b)
rownames(COLDATA)=colnames(b)

sce <- SingleCellExperiment(
    assays = list(
        counts = as.matrix(b),
        logcounts = as.matrix(b)
    ), 
    colData = COLDATA
)

rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
plotPCA(sce)

sce <- sc3(sce, ks = 2:6, biology = TRUE,gene_filter = FALSE)

sc3_interactive(sce)


