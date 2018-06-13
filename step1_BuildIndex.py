################################################
#Author: Feng Zhang
#Email: 15110700005@fudan.edu.cn
#Date: 2018.06.12
################################################


import sys
import cPickle as pickle

print '''
$1 NETWORK_PATH (col1: gene; col2: gene)
$2 EXP_POOL_PATH (row1: header; col 1: gene_name; other_cols: exp value)
$3 OUT_DATA_PATH
$4 PERCENT
'''

NETWORK_PATH=sys.argv[1]
EXP_POOL_PATH=sys.argv[2]
OUT_DATA_PATH=sys.argv[3]
PERCENT=float(sys.argv[4])
#LEN_LIMIT=20

class Person_PCC_Data:
    def __init__(self, EXP, POOL_LENGTH, PCC_POOL,header):
        self.EXP = EXP
        self.POOL_LENGTH = POOL_LENGTH
        self.PCC_POOL = PCC_POOL
        self.header = header


EDGE=set()
POINT=set()
fnet=open(NETWORK_PATH)
for line in fnet:
    seq=line.rstrip().split('\t')
    if line[0]!='#':
        p1=seq[0]
        p2=seq[1]
        edge=[p1,p2]
        edge.sort()
        edge=':'.join(edge)
        EDGE.add(edge)
        POINT.add(p1)
        POINT.add(p2)
fnet.close()

EXP={}
EXP_GENE=set()
POOL_LENGTH=0
fpool=open(EXP_POOL_PATH)
header=fpool.readline().rstrip().split('\t')
for line in fpool:
    #try:
        seq=line.rstrip().split('\t')
        if seq[0] in POINT:
            EXP_GENE.add(seq[0])
            tmp=[]
            for one in seq[1:]:
                try:
                    tmp.append(float(one))
                except Exception as e:
                    print one
                    tmp.append(0.0)

            POOL_LENGTH=len(tmp)
            EXP[seq[0]]=tmp

            if len(header)==len(EXP[seq[0]])+1:
                header=header[1:]


fpool.close()
print POOL_LENGTH



EDGE=set()
POINT=set()
fnet=open(NETWORK_PATH)
for line in fnet:
    seq=line.rstrip().split('\t')
    if line[0]!='#':
        p1=seq[0]
        p2=seq[1]
        if p1 in EXP_GENE and p2 in EXP_GENE:
            edge=[p1,p2]
            edge.sort()
            edge=':'.join(edge)
            EDGE.add(edge)
            POINT.add(p1)
            POINT.add(p2)
fnet.close()

LEN_LIMIT= POOL_LENGTH *  float(PERCENT) #0.2 #20

from scipy import stats
PCC_POOL={}
NETGENE=set()
for edge in EDGE:
    ps=edge.split(':')
    p1=ps[0]
    p2=ps[1]
    #try:
    if EXP[p1].count(0) <= (POOL_LENGTH-LEN_LIMIT) and EXP[p2].count(0) <= (POOL_LENGTH-LEN_LIMIT):
        p1_exp=[]
        p2_exp=[]
        i=0
        while i<POOL_LENGTH:
            exp_p1=EXP[p1][i]
            exp_p2=EXP[p2][i]
            if (exp_p1*exp_p2) !=0:
                p1_exp.append(exp_p1)
                p2_exp.append(exp_p2)
            i+=1
        if len(p1_exp)>=LEN_LIMIT:
            ####################################################
            pcc=stats.pearsonr(p1_exp,p2_exp)[0]
            #pcc=stats.spearmanr(p1_exp,p2_exp)[0]
            ####################################################
            if abs(pcc)!=1:
                PCC_POOL[edge]=[pcc,float(len(p1_exp))]
                NETGENE.add(p1)
                NETGENE.add(p2)
    #except Exception as e:
    #    pass
OUTEXP={}
for gene in EXP:
    if gene in NETGENE:
        OUTEXP[gene]=EXP[gene]

print len(PCC_POOL)
fo=open(OUT_DATA_PATH,'w')
data=Person_PCC_Data( OUTEXP, POOL_LENGTH, PCC_POOL, header)
pickle.dump(data,fo)
fo.close()



