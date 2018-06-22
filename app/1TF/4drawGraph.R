library(stringr)
library(igraph)

top10=read.table('top10.txt',row.names=1,header=T)


cluster_list = unique(top10[,6])


pdf('GRAPH.pdf',width=10,height=10)


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
    
    
    plot(main=as.character(this_cluster), g, vertex.label.cex=2,edge.width=1, vertex.size=5, vertex.label.dist=2, vertex.label.color = "black")
      
    }
dev.off()




