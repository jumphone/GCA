import numpy as np

ACT_FLAG='Activation'
REP_FLAG='Repression'

#TF_TG_FILE='TF_EGFR.info.txt'
TF_TG_FILE='trrust_rawdata.human.tsv'
GCA_OUTPUT='GSE70630'
EDGE_SCORE_CUTOFF=0.01

######################################
######################################

fi=open(TF_TG_FILE)
TFTG={}
TF=set()
TG=set()
for line in fi:
    seq=line.rstrip().split('\t')
    if seq[2]==ACT_FLAG or seq[2]==REP_FLAG :
        this_tf=seq[0]
        this_tg=seq[1]
        TF.add(this_tf)
        TG.add(this_tg)
        if seq[2]==ACT_FLAG:
            this_mode =1
        else:
            this_mode= -1
       
        if this_tf in TFTG:
            TFTG[this_tf][this_tg]=this_mode
        else:
            TFTG[this_tf]={}
            TFTG[this_tf][this_tg]=this_mode 
fi.close()

##############################################

EDGE_SCORE_FILE=GCA_OUTPUT+'/Score_summary.txt'
fi=open(EDGE_SCORE_FILE)
EDGE_SCORE={}
for line in fi:
    seq=line.rstrip().split('\t')
    this_edge_score=float(seq[1])
    if this_edge_score >= EDGE_SCORE_CUTOFF:
        EDGE_SCORE[seq[0]] = this_edge_score  
    
fi.close()

################################################

def SINGLE(this_edge, this_tf, this_tg, this_mode,output_file):
 
    ##############################
    edge_mix_file=GCA_OUTPUT+'/'+this_edge+'.summary.txt'    
    fi=open(edge_mix_file)
    fi.readline()
    info_list=[]
    lambda_list=[]
    for line in fi:
        seq=line.rstrip().split('\t')
        info_list.append([float(seq[1]),seq[4],float(seq[0]) ]) 
        lambda_list.append(float(seq[0]))
    info_list.sort()
    tf_pos_clust='-1'
    if this_mode==1:
        if info_list[-1][2] !=max(lambda_list):
            tf_pos_clust = info_list[-1][1]
        #tf_neg_clust = info_list[0][1]
    elif this_mode==-1:
        if info_list[0][2] !=max(lambda_list):
            tf_pos_clust = info_list[0][1]
        #tf_neg_clust = info_list[-1][1]
    
    fi.close()

    #print tf_pos_clust,tf_neg_clust
    #################################
    edge_cluster_file=GCA_OUTPUT+'/'+this_edge+'.cluster.txt'  
    fi=open(edge_cluster_file) 
    header=fi.readline().rstrip().split('\t')
    tg_index=header.index(this_tg)
    tg_noNA_exp=[]
    for line in fi:
        seq=line.rstrip().split('\t')
        if seq[1]!='NA':
            tg_noNA_exp.append(float(seq[tg_index]))
    tg_cutoff = np.median(tg_noNA_exp)
    #print tg_cutoff
    fi.close()
    
    ###################################
    edge_mode_file = output_file
    fo=open(edge_mode_file,'w')
    fi=open(edge_cluster_file)
    fo.write(fi.readline().rstrip()+'\tMode\n')
    for line in fi:
        seq=line.rstrip().split()
        
        this_tf_flag = 0
        this_tg_flag = 0
        if seq[1] !='NA':
            if float(seq[tg_index]) > tg_cutoff:
                this_tg_flag=1
            elif float(seq[tg_index]) < tg_cutoff:
                this_tg_flag=-1
            if seq[3] == tf_pos_clust:
                this_tf_flag=1
            else: #seq[3] == tf_neg_clust:
                this_tf_flag=-1 
        output_mode='NA'
        #print this_tf_flag, this_tg_flag
        if this_tf_flag* this_tg_flag !=0:
            #print this_tf_flag, this_tg_flag
            if this_tg_flag ==1 and this_tf_flag ==1:
                output_mode='TG_HI;TF_POS;'
            elif this_tg_flag ==1 and this_tf_flag ==-1:
                output_mode='TG_HI;TF_NEG;'
            elif this_tg_flag ==-1 and this_tf_flag ==1:
                output_mode='TG_LW;TF_POS;'
            elif this_tg_flag ==-1 and this_tf_flag ==-1:
                output_mode='TG_LW;TF_NEG;'
        fo.write(line.rstrip()+'\t'+output_mode+'\n')
   
    fi.close()    

########################################

OUTPUT='MODE'
fo=open(OUTPUT+'.summary.txt','w')
for this_edge in EDGE_SCORE:
    p1p2=this_edge.split('.And.')
    p1=p1p2[0]
    p2=p1p2[1]
    
    if p1 in TFTG:
       if p2 in TFTG[p1]:
          this_tf=p1
          this_tg=p2
          this_mode=TFTG[p1][p2]
          if this_mode==1:
              mode_tag='A'
          else:
              mode_tag='R'
          output_file=OUTPUT+'/'+'TF_'+p1+'.TG_'+p2+'.mode.'+str(mode_tag)+'.txt'
          fo.write(OUTPUT+'/'+'TF_'+p1+'.TG_'+p2+'.mode.'+str(mode_tag)+'.txt\n')
          SINGLE(this_edge, this_tf, this_tg, this_mode,output_file)

    elif p2 in TFTG:
       if p1 in TFTG[p2]:
          this_tf=p2
          this_tg=p1
          this_mode=TFTG[p2][p1]
          if this_mode==1:
              mode_tag='A'
          else:
              mode_tag='R'
          output_file=OUTPUT+'/'+'TF_'+p2+'.TG_'+p1+'.mode.'+str(mode_tag)+'.txt'
          fo.write(OUTPUT+'/'+'TF_'+p2+'.TG_'+p1+'.mode.'+str(mode_tag)+'.txt\n')
          SINGLE(this_edge, this_tf, this_tg, this_mode,output_file)


  


