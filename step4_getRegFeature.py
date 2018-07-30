import numpy as np
import subprocess
import sys



###Config#######################################

ACT_FLAG='Activation'
REP_FLAG='Repression'
EDGE_SCORE_CUTOFF=0.01

############################################

print('$1 MODE, $2 GCA_OUTPUT, $3 OUTPUT')
print('')
TF_TG_FILE = sys.argv[1]
GCA_OUTPUT = sys.argv[2]
OUTPUT=sys.argv[3]

###############################

#TF_TG_FILE='trrust_rawdata.human.tsv'
#GCA_OUTPUT='GSE70630'
#OUTPUT='MODE'

OUTPUT_TMP=OUTPUT+'.TMP'
subprocess.Popen('mkdir '+OUTPUT_TMP,shell=True).wait()

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
    tf_index=header.index(this_tf)
    tg_noNA_exp=[]
    tf_noNA_exp=[]
    for line in fi:
        seq=line.rstrip().split('\t')
        if seq[1]!='NA':
            tg_noNA_exp.append(float(seq[tg_index]))
            tf_noNA_exp.append(float(seq[tf_index]))
    tg_cutoff = np.median(tg_noNA_exp)
    tf_cutoff = np.median(tf_noNA_exp)
    #print tg_cutoff
    fi.close()
    
    ###################################
    edge_mode_file = output_file
    fo=open(edge_mode_file,'w')
    fi=open(edge_cluster_file)
    fo.write(fi.readline().rstrip()+'\tMode\n')
    for line in fi:
        seq=line.rstrip().split()
        output_mode='NA'
        
        if seq[1] !='NA':
            if this_mode==1:
                if float(seq[tg_index]) > tg_cutoff and float(seq[tf_index]) > tf_cutoff:
                    output_mode='TG_HI;TF_POS;'
                if float(seq[tg_index]) < tg_cutoff and float(seq[tf_index]) < tf_cutoff:
                    output_mode='TG_HI;TF_POS;'
            if this_mode==-1:
                if float(seq[tg_index]) > tg_cutoff and float(seq[tf_index]) < tf_cutoff:
                    output_mode='TG_HI;TF_POS;'
                if float(seq[tg_index]) < tg_cutoff and float(seq[tf_index]) > tf_cutoff:
                    output_mode='TG_HI;TF_POS;'            

        fo.write(line.rstrip()+'\t'+output_mode+'\n')
   
    fi.close()    

########################################
fo=open(OUTPUT_TMP+'.summary.txt','w')
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
          output_file=OUTPUT_TMP+'/'+'TF_'+p1+'.TG_'+p2+'.mode.'+str(mode_tag)+'.txt'
          fo.write(OUTPUT_TMP+'/'+'TF_'+p1+'.TG_'+p2+'.mode.'+str(mode_tag)+'.txt\n')
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
          output_file=OUTPUT_TMP+'/'+'TF_'+p2+'.TG_'+p1+'.mode.'+str(mode_tag)+'.txt'
          fo.write(OUTPUT_TMP+'/'+'TF_'+p2+'.TG_'+p1+'.mode.'+str(mode_tag)+'.txt\n')
          SINGLE(this_edge, this_tf, this_tg, this_mode,output_file)

        
        
fo.close()        

############################################


fi=open(OUTPUT_TMP+'.summary.txt')
fo=open(OUTPUT,'w')


TFTG={}
for line in fi:
    input_file=line.rstrip()
    tag='TF_'+input_file.split('TF_')[1].split('.txt')[0]

    tf_mode=input_file.split('mode.')[1].split('.txt')[0]
    #print input_file 
    tag=tag.replace('_','.')
    #print tag
    ftmp=open(input_file)
    ftmp.readline()
    out=[tag+'.TG.HI']
    out_lw=[tag+'.TG.LW']
    header=[]
    for l in ftmp:
        sss=l.rstrip().split('\t') 
        header.append(sss[0])
        if sss[-1]=='TG_HI;TF_POS;':
            out.append('1')
        else:
            out.append('0') 

        if sss[-1] == 'TG_LW;TF_POS;':
            out_lw.append('1')
        else:
            out_lw.append('0')
    #if tf_mode=='A':
    TFTG[tag+'.TG.HI']='\t'.join(out)
    #elif tf_mode=='R':
    TFTG[tag+'.TG.LW']='\t'.join(out_lw)
 

fo.write('\t'.join(header)+'\n')
for tag in TFTG:
    fo.write(TFTG[tag]+'\n')


fo.close()


