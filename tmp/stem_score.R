exp_data=read.table('GSE70630_OG_processed_data_v2.txt.cleaned.txt', header=T, row.names=1, check.names=FALSE,sep='\t')
stem_non_EGFR=read.table('STEM.txt',row.names=1)
ST=which(rownames(exp_data) %in% stem_non_EGFR[,1])

s_exp_data=apply(exp_data,1,scale)
s_s_exp_data=apply(s_exp_data,1,scale,scale=F,center=T)

stem_score=apply(s_s_exp_data[ST,], 2, mean)
names(stem_score)=colnames(exp_data)


write.table(stem_score,file='stem_score.txt',row.names=T,col.names=F,sep='\t',quote=F)


