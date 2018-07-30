exp_data=read.table('PDN_Microglia_RawData.txt',header=T,row.names=1)
ident_data=read.table('IDENT.txt')

all_gene=row.names(exp_data)

IDENT=unique(ident_data[,2])
OUT=c()
j=1
while(j<=length(IDENT)){
    this_ident=IDENT[j]
    this_exp_data=exp_data[,which(colnames(exp_data) %in% ident_data[which(ident_data[,2]==this_ident),1] )]
    this_ident_num=length(this_exp_data[1,])
    out=c()
    i=1
    while(i<=length(all_gene)){
        #this_gene=all_gene[i]
        pos_num=length(which(t(this_exp_data[i,])>0))
        out=c(out,pos_num/this_ident_num)
        i=i+1}
    OUT=cbind(OUT,out) 
    print(j)
    j=j+1}

rownames(OUT)=all_gene
colnames(OUT)=as.character(IDENT)

write.table(OUT,file='IDENT_GENE_SIZE.txt',sep='\t',quote=F,row.names=T,col.names=T)






