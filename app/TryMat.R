exp_data=read.table('../GSE70630_OG_processed_data_v2.txt.cleaned.txt',header=T,row.names=1)
exp_data[is.na(exp_data)]=0

value=as.numeric(as.matrix(exp_data))
value=value[which(value>0)]
length(value)

gene_num=length(exp_data[,1])
cell_num=length(exp_data[1,])

SEED=1

while(SEED<=10){
set.seed(SEED)
exp_value=c()

i=1
while(i<=cell_num){
    exp_value=exp_data[,i]
    exp_value_pos = exp_value[which(exp_value>0)]
    #m=mean(exp_value_pos)
    #s=sd(exp_value_pos)
    #r_exp_value_pos=rnorm(length(exp_value_pos),mean=m,sd=d)
    r_exp_value_pos=sample(value,size=length(exp_value_pos))
    #exp_data[,i]=rep(0,gene_num)
    #exp_data[sample(c(1:gene_num),size=length(exp_value_pos)),i]=r_exp_value_pos
    exp_data[which(exp_value>0),i]=r_exp_value_pos
    i=i+1
    print(i)
    }


write.table(exp_data, file=paste0('random_exp_data_',as.character(SEED),'.txt'),row.names=T,col.names=T,sep='\t',quote=F)
SEED=SEED+1
}
