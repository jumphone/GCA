################################################
#Author: Feng Zhang
#Email: 15110700005@fudan.edu.cn
#Date: 2018.06.12
################################################

########Load packages############
library(stringr)
library(pastecs)
library(mixtools)
library(parallel)
library(igraph)

########Args############
print('$1 EXP, $2 ZMAT, $3 OUT, $4 CPU, $5 SEED')
args = commandArgs(trailingOnly=TRUE)
REF=args[1]
INPUT=args[2]
OUTPUT=args[3]
CPU=as.numeric(args[4])
RANDOM_SEED=as.numeric(args[5])
TMP_DIR=paste0(OUTPUT)
system(paste0('mkdir ',TMP_DIR))
names(args)=c('EXP','ZMAT','OUT','CPU','SEED')
write.table(args, file=paste0(TMP_DIR,'/Arguments.txt'),sep='\t',quote=F,row.names=T, col.names=F )
########################


#########Read Data############
exp_data=read.table(REF, header=T, row.names=1, check.names=FALSE,sep='\t')
gene_name=rownames(exp_data)
input_data=read.table(INPUT,header=T,row.names=1,check.names=FALSE,sep='\t')
ROW_NUM=length(input_data[,1])
COL_NUM=length(input_data[1,])
ROW_LABEL=rownames(input_data)
COL_LABEL=colnames(input_data)
########################

######Single Thread Funtion##################
SINGLE = function(i){
    print(i)
    set.seed(RANDOM_SEED)
    this_row_label=ROW_LABEL[i]
    tmp=t(input_data[i,])
    
    ###########Remove NA and Outliers#############
    tmp=tmp[which(!is.na(tmp))]
    IQR=quantile(tmp,0.75)-quantile(tmp,0.25)
    UP=quantile(tmp,0.75)+1.5*IQR
    DW=quantile(tmp,0.25)-1.5*IQR
    tmp=tmp[which(tmp<UP & tmp >DW)]
    ############################################
    ###########Find peak#############
    D=density(tmp)
    PEAK_PIT=extract(turnpoints(D$y),length(D$y),peak=1,pit=-1)
    MEAN=D$x[which(PEAK_PIT==1)]
    PEAK_NUM=length(which(PEAK_PIT==1))

    tmp_out_path=paste0(TMP_DIR,'/',this_row_label)
    ori_data = t(input_data[i,])
    ############################################
    if(PEAK_NUM>1){       
        ###########Mixtools & assign cluster#############
        mix1=normalmixEM(tmp,mu=MEAN,mean.constr=MEAN,maxit=10000)        
        clust_out = rep(0,length(ori_data))
        j=1
        while(j<=length(ori_data)){
            this_z = ori_data[j]
            if(!is.na(this_z) & this_z < UP & this_z > DW){
                this_d = dnorm(this_z, mean=mix1$mu, sd=mix1$sigma )*mix1$lambda
                this_c=1
                while(this_c <=length(this_d)){
                    if(this_d[this_c] >=max(this_d)){clust_out[j]=this_c}
                    this_c=this_c+1
                    }                
                }
            j=j+1
            }     
        ############################################
        #if(length(unique(clust_out[which(clust_out!=0)]))>1  ) { 
        if(1==1){ 
            ###########Draw figures#############
            pdf(paste0(tmp_out_path,'.pdf'),width=10,height=10)
            par(mfrow=c(2,2))
            plot(D,main=this_row_label)
            abline(v=MEAN,col='red',lty=3)
            plot.mixEM(mix1,whichplots=2,breaks=50)
            COL=clust_out+1

            ############################  
            all_cell_num = length(COL_LABEL)
            this_tag_cell_num = length(tmp)
            this_pie_data=c(1-this_tag_cell_num/all_cell_num)
            this_pie_col=c(0)
            mix_index=1
            while(mix_index<=length(mix1$lambda)){
                this_pie_data=c(this_pie_data, mix1$lambda[mix_index] * this_tag_cell_num/all_cell_num )
                this_pie_col=c(this_pie_col,mix_index)
                mix_index=mix_index+1}       
            this_pie_col=this_pie_col+1
            this_pie_col=palette()[this_pie_col]
            this_pie_col[1]='grey80'
            pie(this_pie_data, labels=as.character(round(this_pie_data,2)),col=this_pie_col, radius = 0.9, main='Proportion Estimation (All)')
            ##############################        
            pie(mix1$lambda, labels=as.character(round(mix1$lambda,2)),col=(c(1:length(mix1$lambda))+1), radius = 0.9, main='Proportion Estimation (nonNA)')
            ################################
            
            cell_index=c(1:length(ori_data))
            this_pch=rep(16,length(ori_data)) 
 
            this_pch[which(COL==1)]=3
          
            plot(ori_data, cell_index, col=COL, pch=this_pch, xlab='z_value', main='All')  #,xlim=c(DW,UP))
            abline(v=UP,col='black',lty=3)
            abline(v=DW,col='black',lty=3)

            plot(ori_data, cell_index, col=COL, pch=this_pch  ,xlim=c(DW,UP),xlab='z_value', main='No outlier')
            abline(v=UP,col='black',lty=3)
            abline(v=DW,col='black',lty=3)
            
         
            p1p2=unlist(strsplit(this_row_label, ".And."))
            p1=p1p2[1]
            p2=p1p2[2]
            p1_exp=t(exp_data[which(gene_name == p1),])
            p2_exp=t(exp_data[which(gene_name == p2),])

            this_xlim = c(min(p1_exp),max(p1_exp))
            this_ylim = c(min(p2_exp),max(p2_exp))
            
            this_v = which(p1_exp!=0 & p2_exp!=0)  
            this_v_out = which( !(ori_data < UP & ori_data > DW ) )
            this_col=rep('black',length(p1_exp))
            this_col[this_v_out]='grey'
            this_pch=rep(16,length(p1_exp))
            this_pch[this_v_out]=3
            this_pcc=round(cor(p1_exp[this_v], p2_exp[this_v]),2)
            plot(p1_exp[this_v], p2_exp[this_v], xlab=p1, ylab=p2, xlim=this_xlim, ylim=this_ylim,main=paste0('All, N=',as.character(length(this_v)),', PCC=',as.character(this_pcc)),pch = this_pch[this_v] , col=this_col[this_v])
            

            ############################
            plot(p1_exp[this_v], p2_exp[this_v], xlab=p1, ylab=p2, xlim=this_xlim, ylim=this_ylim, main='With color', pch = this_pch[this_v] , col=this_col[this_v])
            for(this_cluster_index in unique(clust_out[which(clust_out!=0)])){
                tmp_cell_index = which(clust_out == this_cluster_index)
                
                if(length(tmp_cell_index) >= 1 ){
                    tmp_p1_exp=p1_exp[tmp_cell_index]
                    tmp_p2_exp=p2_exp[tmp_cell_index]
                    par(new=TRUE)
                    plot(tmp_p1_exp, tmp_p2_exp, xlim=this_xlim, ylim=this_ylim,main='',xlab='',ylab='',pch=16,col=this_cluster_index+1)
                    }
                }
            ############################
            this_v_nout = which( (ori_data < UP & ori_data > DW ) )
            this_nout_xlim = c(min(p1_exp[this_v_nout ]),max(p1_exp[this_v_nout ]))
            this_nout_ylim = c(min(p2_exp[this_v_nout ]),max(p2_exp[this_v_nout ]))

            for(this_cluster_index in unique(clust_out[which(clust_out!=0)])){
                tmp_cell_index = which(clust_out == this_cluster_index)
                
                if(length(tmp_cell_index) >= 1 ){
                    tmp_p1_exp=p1_exp[tmp_cell_index]
                    tmp_p2_exp=p2_exp[tmp_cell_index]
                    this_pcc=round(cor(c(p1_exp[this_v],tmp_p1_exp), c(p2_exp[this_v],tmp_p2_exp)),2)
                    plot(tmp_p1_exp, tmp_p2_exp, xlab=p1, ylab=p2, xlim=this_nout_xlim, ylim=this_nout_ylim,main=paste0('Cluster',as.character(this_cluster_index),', N=',as.character(length(tmp_p1_exp)),', AddPCC=',as.character(this_pcc)),pch=16,col=this_cluster_index+1)                   
                    }  
                }


            dev.off()
            ############################
            ###########Write files#################
            mix1.cluster=c(1:length(mix1$lambda))
            mix1.color=palette()[c(1:length(mix1$lambda))+1]
            tmp_out=cbind(mix1$lambda,mix1$mu,mix1$sigma,mix1.color,mix1.cluster)
            colnames(tmp_out)=c('lambda','mu','sigma','color','cluster')         
            write.table(as.matrix(tmp_out),file=paste0(tmp_out_path,'.summary.txt'),sep='\t',quote=F,row.names=F,col.names=T)
            
            clust_out=cbind(COL_LABEL, ori_data, palette()[clust_out+1], clust_out, p1_exp, p2_exp)
            colnames(clust_out)=c('Cell_name','Zvalue','Color','Cluster',p1,p2) 
            write.table(as.matrix(clust_out),file=paste0(tmp_out_path,'.cluster.txt'),sep='\t',quote=F,row.names=F,col.names=T)          
            ############################
            if(length(mix1$lambda)>=2){
            #OUT=c(this_row_label, sort(mix1$lambda,decreasing=T)[2])
                OUT=c(this_row_label, length(tmp) * sort(mix1$lambda / length(COL_LABEL),decreasing=T)[2])
                }
            else{OUT=c(this_row_label,0)}
            return(OUT)
            }
        }
    else{
        if(1==1){
            #######draw#############
            pdf(paste0(tmp_out_path,'.pdf'),width=10,height=10)
            par(mfrow=c(2,2))
            plot(D,main=this_row_label)
            abline(v=MEAN,col='red',lty=3)
            hist(tmp,breaks=50)
            ############################
            cell_index=c(1:length(ori_data))
            this_pch=rep(16,length(ori_data)) 
            this_v_out = which( !(ori_data < UP & ori_data > DW ) )
            this_pch[this_v_out]=3  
            plot(ori_data, cell_index,  pch=this_pch, xlab='z_value', main='All')  #,xlim=c(DW,UP))
            abline(v=UP,col='black',lty=3)
            abline(v=DW,col='black',lty=3)
            plot(ori_data, cell_index, pch=this_pch  ,xlim=c(DW,UP),xlab='z_value', main='No outlier')
            abline(v=UP,col='black',lty=3)
            abline(v=DW,col='black',lty=3)
            ############################
            p1p2=unlist(strsplit(this_row_label, ".And."))
            p1=p1p2[1]
            p2=p1p2[2]
            p1_exp=t(exp_data[which(gene_name == p1),])
            p2_exp=t(exp_data[which(gene_name == p2),])
             
            this_xlim = c(min(p1_exp),max(p1_exp))
            this_ylim = c(min(p2_exp),max(p2_exp))
            
            this_v = which(p1_exp!=0 & p2_exp!=0)  
            this_v_out = which( !(ori_data < UP & ori_data > DW ) )
            this_col=rep('black',length(p1_exp))
            this_col[this_v_out]='grey'
            this_pch=rep(16,length(p1_exp))
            this_pch[this_v_out]=3
            this_pcc=round(cor(p1_exp[this_v], p2_exp[this_v]),2)
            plot(p1_exp[this_v], p2_exp[this_v], xlab=p1, ylab=p2, xlim=this_xlim, ylim=this_ylim,main=paste0('All, N=',as.character(length(this_v)),', PCC=',as.character(this_pcc)),pch = this_pch[this_v] , col=this_col[this_v])

            dev.off()
            ###########Write files#################
            #mix1.cluster=c(1:length(mix1$lambda))
            #mix1.color=palette()[c(1:length(mix1$lambda))+1]
            #tmp_out=cbind(mix1$lambda,mix1$mu,mix1$sigma,mix1.color,mix1.cluster)
            #colnames(tmp_out)=c('lambda','mu','sigma','color','cluster')    
            tmp_out='None'     
            write.table(as.matrix(tmp_out),file=paste0(tmp_out_path,'.summary.txt'),sep='\t',quote=F,row.names=F,col.names=T)
            
            clust_out=cbind(COL_LABEL, ori_data, p1_exp, p2_exp)
            colnames(clust_out)=c('Cell_name','Zvalue',p1,p2) 
            write.table(as.matrix(clust_out),file=paste0(tmp_out_path,'.cluster.txt'),sep='\t',quote=F,row.names=F,col.names=T)          
            ############################
            #OUT=c(this_row_label, sort(mix1$lambda,decreasing=T)[2])
            OUT=c(this_row_label, 0)
            return(OUT)            
            }
        }  
    }
#######################################

#######################################
#######################################
#######################################
RUN = mclapply(1:ROW_NUM, SINGLE, mc.cores=CPU)
#######################################
#######################################
#######################################


SECOND_LAMBDA=c('Tag','Second_Lambda')
for(one in RUN){
    if(length(one)>1){
        SECOND_LAMBDA=cbind(SECOND_LAMBDA,one)
        }
    }


SECOND_LAMBDA = t(SECOND_LAMBDA )
OUT_SECOND_LAMBDA = as.numeric(SECOND_LAMBDA[c(2:length(SECOND_LAMBDA[,1])),2])
OUT_SECOND_LAMBDA=as.matrix(OUT_SECOND_LAMBDA)
rownames(OUT_SECOND_LAMBDA) = SECOND_LAMBDA[c(2:length(SECOND_LAMBDA[,1])),1]

O=order( OUT_SECOND_LAMBDA[,1], decreasing=T)
OUT_SECOND_LAMBDA=as.matrix(OUT_SECOND_LAMBDA[O,])
#OUT_SECOND_LAMBDA=as.matrix(round(OUT_SECOND_LAMBDA,2))
SC_OUT_SECOND_LAMBDA=as.matrix(format(OUT_SECOND_LAMBDA,digits = 2, scientific = TRUE))
write.table(SC_OUT_SECOND_LAMBDA,file=paste0(TMP_DIR,'/Score_summary.txt'),sep='\t',quote=F,row.names=T,col.names=F )


####GENE_RANK##################
i=1
RANK_GENE=c()
while(i <= length(OUT_SECOND_LAMBDA[,1])){
    this_tag = rownames(OUT_SECOND_LAMBDA)[i]    
    this_second_lambda = SECOND_LAMBDA[i,1]
    p1p2=unlist(strsplit(this_tag, ".And."))
    p1=p1p2[1]
    p2=p1p2[2]
    RANK_GENE=c(RANK_GENE,p1,p2)
    i=i+1} 
RANK_GENE=unique(RANK_GENE)

RANK_GENE_KSP=c()
RANK_GENE_KSP_NUM=c()

SCORE_LIST=cbind(OUT_SECOND_LAMBDA[,1],OUT_SECOND_LAMBDA[,1],OUT_SECOND_LAMBDA[,1])
j=1
while(j<=length(OUT_SECOND_LAMBDA[,1])){
    this_tag = rownames(OUT_SECOND_LAMBDA)[j]    
    this_second_lambda = OUT_SECOND_LAMBDA[j,1]
    p1p2=unlist(strsplit(this_tag, ".And."))
    p1=p1p2[1]
    p2=p1p2[2]
    SCORE_LIST[j,2]=p1
    SCORE_LIST[j,3]=p2
    j=j+1
    }

i=1
while(i<=length(RANK_GENE)){
    print(i)
    this_rank_gene=RANK_GENE[i]
    this_rank_gene_second_lambda=as.numeric(SCORE_LIST[which(SCORE_LIST[,2]==this_rank_gene | SCORE_LIST[,3]==this_rank_gene),1])
    this_rank_gene_background_second_lambda=as.numeric(SCORE_LIST[which(SCORE_LIST[,2]!=this_rank_gene & SCORE_LIST[,3]!=this_rank_gene),1])
    this_rank_gene_p = ks.test(this_rank_gene_background_second_lambda, this_rank_gene_second_lambda,alternative='greater')$p.value
    RANK_GENE_KSP=c(RANK_GENE_KSP,this_rank_gene_p)
    RANK_GENE_KSP_NUM=c(RANK_GENE_KSP_NUM,length(this_rank_gene_second_lambda))
    i=i+1
    }

names(RANK_GENE_KSP)=RANK_GENE
RANK_GENE_KSP_NUM = RANK_GENE_KSP_NUM[order(RANK_GENE_KSP)]
RANK_GENE_KSP = RANK_GENE_KSP[order(RANK_GENE_KSP)]
RANK_GENE_KSP_OUT=cbind(RANK_GENE_KSP_NUM,RANK_GENE_KSP)
rownames(RANK_GENE_KSP_OUT)=names(RANK_GENE_KSP)
write.table(RANK_GENE_KSP_OUT, file=paste0(TMP_DIR,'/Pvalue_summary.txt'),sep='\t',quote=F,row.names=T, col.names=F )

#######################

##########Draw graph#####################
set.seed(RANDOM_SEED)
NET = cbind(rep('tag',length(OUT_SECOND_LAMBDA)),rep('tag',length(OUT_SECOND_LAMBDA)))
i=1
while(i<=length(OUT_SECOND_LAMBDA[,1])){
    p1p2 = unlist(strsplit(rownames(OUT_SECOND_LAMBDA)[i], ".And."))
    p1 = p1p2[1]
    p2 = p1p2[2]
    NET[i,1]=p1
    NET[i,2]=p2
    i=i+1}
g <- make_graph(t(NET),directed = FALSE)
colors <- colorRampPalette(c('white','lightpink','indianred1',"red1", "red2", "red3", "red4",'darkred','darkred','darkred','darkred'))(51)
E(g)$color = colors[as.integer(OUT_SECOND_LAMBDA[,1] * 100)+1]

node.size=setNames( (1-RANK_GENE_KSP)*3,names(RANK_GENE_KSP))

pdf(paste0(TMP_DIR,'/G.pdf'),width=20,height=20)
plot( c(1:51)/100,c(1:51)/100, col=colors, ylab='Score', xlab='Score', pch=16,cex=5,lwd=5,type='p',main='Edge Color Key')

#l <- layout_with_fr(g)
#plot(main='All', g, layout=layout_with_fr, vertex.label.cex=0.5, vertex.size=as.matrix(node.size), vertex.label.dist=0, vertex.label.color = "black",vertex.frame.color = "white",vertex.color = "gold2")
plot(main='All', g, vertex.label.cex=0.5, vertex.size=as.matrix(node.size), vertex.label.dist=0, vertex.label.color = "black",vertex.frame.color = "white",vertex.color = "gold2")

V(g)$comp <- components(g)$membership
i=1
while(i<=max(V(g)$comp)){
    this_subg = induced_subgraph(g,V(g)$comp==i)
    #l <- layout_with_fr(this_subg)
    #plot(main=paste0('SubGraph',as.character(i)),this_subg, layout=layout_with_fr, vertex.label.cex=0.5, vertex.size=as.matrix(node.size), vertex.label.dist=0, vertex.label.color = "black",vertex.frame.color = "white",vertex.color = "gold2")
    if(length(E(this_subg))>=2){
        plot(main=paste0('SubGraph',as.character(i)),this_subg, vertex.label.cex=0.5, vertex.size=as.matrix(node.size), vertex.label.dist=0, vertex.label.color = "black",vertex.frame.color = "white",vertex.color = "gold2")
        }
    i=i+1}
dev.off()
#######################

WARNINGS=warnings();

#########Output HTML#############
OUT_HTML=c('<center><header><h1>Result</h1></header>')
OUT_HTML=c(OUT_HTML, '<a href="G.pdf"> Graph </a> &nbsp&nbsp&nbsp <a href="Arguments.txt">Args</a></br></br>')
#######################
OUT_HTML=c(OUT_HTML,'<table width="80%" align="center"><tr><td width="40%" valign="top">')
OUT_HTML=c(OUT_HTML,'<table border="1", width="100%"><tr><td>NO.</td><td>Tag</td><td>Num</td><td>Pvalue</td></tr>')
i=1
while(i<=length(RANK_GENE_KSP)){
    this_tag = names(RANK_GENE_KSP)[i]
    this_pvalue = format(RANK_GENE_KSP[i],digits = 2, scientific = TRUE)
    this_num = RANK_GENE_KSP_NUM[i]
    this_out = paste0('<tr><td>',as.character(i),'</td><td>',this_tag,'</td><td>',as.character(this_num), '</td><td>',as.character(this_pvalue),'</td></tr>')
    OUT_HTML=c(OUT_HTML,this_out)
    i=i+1
    }
OUT_HTML=c(OUT_HTML,'<td><a href="Pvalue_summary.txt">TXT</a></td></table>')

OUT_HTML=c(OUT_HTML,'</td>' )
######################
OUT_HTML=c(OUT_HTML,'<td width="5%" valign="top"></td>')
OUT_HTML=c(OUT_HTML,'<td width="55%" valign="top">')
OUT_HTML=c(OUT_HTML,'<table border="1", width="100%"><tr><td>NO.</td><td>Tag</td><td>Score</td><td>Figure</td><td>MixInfo</td><td>Cluster</td></tr>')
i=1
while(i <=length(SC_OUT_SECOND_LAMBDA[,1])){
    this_tag = rownames(SC_OUT_SECOND_LAMBDA)[i]
    second_lambda = SC_OUT_SECOND_LAMBDA[i,1]
    this_out = paste0('<tr><td>',as.character(i),'</td><td>',this_tag,'</td><td>', as.character(second_lambda),'</td><td>','<a href="',this_tag,'.pdf"> pdf </a>','</td><td>' ,'<a href="',this_tag,'.summary.txt"> txt </a>','</td><td>',' <a href="',this_tag,'.cluster.txt"> txt </a></td> </tr>')
    OUT_HTML=c(OUT_HTML,this_out)
    i=i+1
    }
OUT_HTML=c(OUT_HTML,'<td><a href="Score_summary.txt">TXT</a></td></table>')
OUT_HTML=c(OUT_HTML,'</td></tr></table>' )

############################
OUT_HTML=c(OUT_HTML,'</center>')
write.table(OUT_HTML,file=paste0(TMP_DIR,'/index.html'),sep='\t',quote=F,row.names=F,col.names=F )
############################



