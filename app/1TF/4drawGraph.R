library(stringr)
library(igraph)

top10=read.table('top10.txt',row.names=1,header=T)


cluster_list = unique(top10[,6])


pdf('GRAPH.pdf',width=10,height=10)


for(this_cluster in cluster_list){
    
    this_cluster_info=top10[which(top10[,6]==this_cluster),]
    
    NET = cbind(rep('tag',length(this_cluster_info[,1])),rep('tag',length(this_cluster_info[,1])))  
    EDGE_COLOR = c()
    i=1
    while(i<=length(this_cluster_info[,1])){
        
        this_tftg_info = unlist(strsplit(as.character(this_cluster_info[i,7]), '.',fixed=TRUE))   
        this_tf=this_tftg_info[2]
        this_tg=this_tftg_info[4]
        this_mode=this_tftg_info[6]
        this_tg_exp=this_tftg_info[8]
        if(this_mode=='A'){EDGE_COLOR = c(EDGE_COLOR ,'red')}
        else{EDGE_COLOR = c(EDGE_COLOR ,'blue')}
        print(this_tftg_info)
        
        NET[i,1]=this_tf
        NET[i,2]=this_tg
        
        
        i=i+1}
        g <- make_graph(t(NET),directed = TRUE)
        
        
        E(g)$color = EDGE_COLOR
        
        plot(main='', g, vertex.label.cex=1,edge.width=1, vertex.size=2, vertex.label.dist=0, vertex.label.color = "black",vertex.frame.color = "white",vertex.color = "gold2")
      
    }
dev.off()




