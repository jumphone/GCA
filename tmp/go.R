library(gplots)
a=as.matrix(read.table('MODE_MAT.txt',header=T,row.names=1))
SUM=apply(a,2,sum)
B= which(SUM>=3)
b=a[,B]
RSUM=apply(b,1,sum)
b=b[which(RSUM>0),]
heatmap(b,Rowv=T,Colv=F,scale='row',labCol='',margins=c(10,10))


library(Seurat)

