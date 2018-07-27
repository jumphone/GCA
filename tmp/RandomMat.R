exp_data=read.table('GSE70630_OG_processed_data_v2.txt.cleaned.txt',header=T,row.names=1)
exp_data[is.na(exp_data)]=0

gene_num=length(exp_data[,1])
cell_num=length(exp_data[1,])

exp_value=c()

i=1
while(i<=cell_num){
    exp_value=exp_data[,i]
    exp_value_pos = exp_value[which(exp_value>0)] 
    r_exp_value_pos=sample(exp_value_pos,size=length(exp_value_pos))
    exp_data[which(exp_value>0),i]=r_exp_value_pos
    i=i+1
    print(i)
    }
  
  
write.table(exp_data, file='random_exp_data.txt',row.names=T,col.names=T,sep='\t',quote=F)
