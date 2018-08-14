library(gplots)
a=as.matrix(read.table('MODE_MAT.txt',header=T,row.names=1))
stem_score=read.table('stem_score.txt')[,2]
ac_score=read.table('AC_score.txt')[,2]
oc_score=read.table('OC_score.txt')[,2]
IDENT=read.table('IDENT.txt',row.names=T,header=F,sep='\t')

