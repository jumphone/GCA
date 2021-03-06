
################################################################
# Cluster
################################################################


library(Seurat)

a = as.matrix(read.table('MODE_MAT.txt',header=T,row.names=1))
SUM = apply(a,2,sum)
B = which(SUM>=3)
b = a[,B]
RSUM=apply(b,1,sum)
b=b[which(RSUM>0),]

all_gene=rownames(b)

EXP = CreateSeuratObject(raw.data = b, min.cells = 0, min.genes=0)
EXP=NormalizeData(object = EXP, normalization.method = "LogNormalize", scale.factor = 10000)
EXP = ScaleData(object = EXP,, genes.use = all_gene)

PCNUM=40
EXP <- RunPCA(object = EXP, pc.genes = all_gene, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )
PCElbowPlot(object = EXP,num.pc=PCNUM)

#PCAPlot(object = EXP, dim.1 = 1, dim.2 = 2)

PCUSE=1:4
EXP = RunTSNE(object = EXP, dims.use = PCUSE, do.fast = TRUE,check_duplicates = FALSE )

RES=0.5
EXP <- FindClusters(object = EXP, reduction.type = "pca", dims.use = PCUSE,  resolution = RES, print.output = 0, save.SNN = TRUE,force.recalc =T,algorithm=1,k.param = 10,n.iter=10)
TSNEPlot(object = EXP,do.label=T)

EXP@scale.data=as.matrix(EXP@data)
pbmc=EXP
library(dplyr)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.05, thresh.use = 0.1)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(30, avg_logFC)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,col.low = "grey90", col.mid = "grey75", col.high = "red",cex.row=6 )

pdf('TSNE.pdf',width=5,height=5)
TSNEPlot(object = EXP,do.label=T)
dev.off()
pdf('HEAT.pdf',width=15,height=15)
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE,col.low = "grey90", col.mid = "grey90", col.high = "red",cex.row=6 )
dev.off()

write.table(file='IDENT.txt',EXP@ident,row.names=T,col.names=F,sep='\t',quote=F)
write.table(pbmc.markers ,file='markers.txt',row.names=T,col.names=T,quote=F,sep='\t')
write.table(top10,file='top10.txt',row.names=T,col.names=T,quote=F,sep='\t')


################################################################
#Survival
################################################################

library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)

a=read.table('metadata_tag.txt.GBM',sep='\t')
b=read.table('IDENT.txt')
a=a[which(a[,23] %in% b[,1]),]
b=b[which(b[,1] %in% a[,23]),]
dim(a)
dim(b)

SUR=as.numeric(as.character(a[,19]))
data=cbind(b[,2],SUR)
data[,2]=data[,2]/30
#boxplot(data[,2]~data[,1])

summary

status=rep(1,length(data[,1]))
age=-as.numeric(as.character(a[,18]))
data=cbind(data,status,age)
colnames(data)=c('cluster','time','status','age')
data=data[which(!is.na(data[,2])),]
data=as.data.frame(data)
data=base::unique(data)

km=with(data,Surv(time,status))
km_fit = survfit(Surv(time,status)~1,data=data)
km_cluster_fit = survfit(Surv(time,status) ~ cluster, data=data)
tmp_km_cluster_fit = km_cluster_fit
tmp_km_cluster_fit$lower=tmp_km_cluster_fit$surv
tmp_km_cluster_fit$upper=tmp_km_cluster_fit$surv
pdf('SUR.pdf')
km_cluster_fit
sdf=survdiff(formula = Surv(time, status)~cluster , data =data, rho = 0)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
p.val=format(p.val, scientific = TRUE) 
autoplot(km_cluster_fit,main=paste0('chisq_test_p = ', as.character(p.val)) )
autoplot(tmp_km_cluster_fit,main=paste0('chisq_test_p = ', as.character(p.val)) )
dev.off()

write.table(file='SUR_N.txt',sdf$n,row.names=T,col.names=T,quote=F,sep='\t')

MEAN_SUR=aggregate(time~cluster,data=data,mean)
write.table(file='SUR_MEAN.txt',MEAN_SUR,row.names=T,col.names=T,quote=F,sep='\t')


################################################################
# Draw Graph
################################################################

library(stringr)
library(igraph)

top10=read.table('top10.txt',row.names=1,header=T)

cluster_list = unique(top10[,6])

pdf('GRAPH.pdf',width=15,height=15)
#par(mfrow=c(2,2))
for(this_cluster in cluster_list){   
    this_cluster_info=top10[which(top10[,6]==this_cluster),]    
    NET = cbind(rep('tag',length(this_cluster_info[,1])),rep('tag',length(this_cluster_info[,1])))  
    EDGE_COLOR = c()
    NODE_COLOR = c()
    i=1
    while(i<=length(this_cluster_info[,1])){
        
        this_tftg_info = unlist(strsplit(as.character(this_cluster_info[i,7]), '.',fixed=TRUE))   
        this_tf=this_tftg_info[2]
        this_tg=this_tftg_info[4]
        this_mode=this_tftg_info[6]
        this_tg_exp=this_tftg_info[8]
        if(this_mode=='A'){EDGE_COLOR = c(EDGE_COLOR ,'red')}
        else{EDGE_COLOR = c(EDGE_COLOR ,'blue')}
                
        if(this_tg_exp=='HI'){NODE_COLOR = cbind(NODE_COLOR ,c(this_tg,'red'))}
        else{NODE_COLOR = cbind(NODE_COLOR ,c(this_tg,'blue'))}
        
        NODE_COLOR=as.matrix(NODE_COLOR)
        
        if(! this_tf %in% t(NODE_COLOR)[,1] ){ 
            if(this_tg_exp=='HI' & this_mode=='A'){ NODE_COLOR=cbind(NODE_COLOR ,c(this_tf, 'red'))}
            if(this_tg_exp=='LW' & this_mode=='A'){ NODE_COLOR=cbind(NODE_COLOR ,c(this_tf, 'blue'))}
            if(this_tg_exp=='HI' & this_mode=='R'){ NODE_COLOR=cbind(NODE_COLOR ,c(this_tf, 'blue'))}
            if(this_tg_exp=='LW' & this_mode=='R'){ NODE_COLOR=cbind(NODE_COLOR ,c(this_tf, 'red'))}
        }
        
        print(this_tftg_info)
        
        NET[i,1]=this_tf
        NET[i,2]=this_tg        
        i=i+1}
    
    g <- make_graph(t(NET),directed = TRUE)    
    E(g)$color = EDGE_COLOR
    NODE_COLOR=t(as.matrix(NODE_COLOR))
    #node.color=setNames( t(NODE_COLOR)[,2],t(NODE_COLOR)[,1])
    
    NEW_NODE_COLOR=c()
    VG=as.character(V(g))
    for(vg in names(V(g))){     
        this_node_color= NODE_COLOR[which(NODE_COLOR[,1]==as.character(vg)),2]  
        V(g)[vg]$color<-this_node_color
        }   
    plot(main=as.character(this_cluster), g, vertex.label.cex=1.5,edge.width=1, vertex.size=10, vertex.label.dist=1.3, vertex.label.color = "black")
      
    }
dev.off()



