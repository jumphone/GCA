exp_data=read.table('GSE70630_OG_processed_data_v2.txt.cleaned.txt', header=T, row.names=1, check.names=FALSE,sep='\t')

s_exp_data=apply(exp_data,1,scale)
s_s_exp_data=apply(s_exp_data,1,scale,scale=F,center=T)

stem=read.table('STEM.txt')
ST=which(rownames(exp_data) %in% stem[,1])
stem_score=apply(s_s_exp_data[ST,], 2, mean)
names(stem_score)=colnames(exp_data)
write.table(stem_score,file='stem_score.txt',row.names=T,col.names=F,sep='\t',quote=F)


stem=read.table('AC.txt')
ST=which(rownames(exp_data) %in% stem[,1])
stem_score=apply(s_s_exp_data[ST,], 2, mean)
names(stem_score)=colnames(exp_data)
write.table(stem_score,file='AC_score.txt',row.names=T,col.names=F,sep='\t',quote=F)


stem=read.table('OC.txt')
ST=which(rownames(exp_data) %in% stem[,1])
stem_score=apply(s_s_exp_data[ST,], 2, mean)
names(stem_score)=colnames(exp_data)
write.table(stem_score,file='OC_score.txt',row.names=T,col.names=F,sep='\t',quote=F)



