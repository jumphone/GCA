################################################
#Author: Feng Zhang
#Email: 15110700005@fudan.edu.cn
#Date: 2018.06.12
################################################

import cPickle as pickle
import sys
from scipy import stats
import random,subprocess,os
import multiprocessing

print '''
$1: Index Data
$2: Output file
$3: CPU
'''

INDEX_DATA_FILE=sys.argv[1]
#EXP_MATRIX=sys.argv[2]
OUT_FILE=sys.argv[2]
PROC_LIMIT=int(sys.argv[3])



class Person_PCC_Data:
    def __init__(self,  EXP , POOL_LENGTH , PCC_POOL, header ):
        self.EXP = EXP
        self.POOL_LENGTH = POOL_LENGTH
        self.PCC_POOL = PCC_POOL
        self.header = header

print 'loading...'
fdata = open(INDEX_DATA_FILE)
data = pickle.load(fdata)
fdata.close()
print "loading done !"

header=data.header
NETGENE=set()
for edge in data.PCC_POOL:
    p1=edge.split(':')[0]
    p2=edge.split(':')[1]
    NETGENE.add(p1)
    NETGENE.add(p2)
    

OUTDIR=OUT_FILE+'_TMP_'+str(random.random())
subprocess.Popen('mkdir '+OUTDIR,shell=True).wait()

def SINGLE(p1, p2,p1_old_exp,p2_old_exp, p1_this_list, p2_this_list, pcc_old,pcc_length):
    outfile=OUTDIR+'/'+p1+'_'+p2+'.ssn'
    if os.path.exists(outfile)!=True:
        Z=[]
        j=0
        while j<len(header):
            p1_this=p1_this_list[j]
            p2_this=p2_this_list[j]
            if p1_this * p2_this>0: #and data.PCC_POOL[edge][1]>=5:
                p1_new_exp=p1_old_exp+[p1_this]
                p2_new_exp=p2_old_exp+[p2_this]
                pcc_new = stats.pearsonr(p1_new_exp,p2_new_exp)[0]
                #pcc = data.PCC_POOL[edge][0]
                delta_pcc = pcc_new - pcc_old
                #z=delta_pcc /( (1-pcc**2)/float(data.PCC_POOL[edge][1]-1) )
                z=delta_pcc /( (1-pcc_old**2)/float(pcc_length -1) )
            else:
                z='NA'
            Z.append(str(z))
            j+=1
        fo=open(outfile,'w')
        fo.write(p1+'.And.'+p2+'\t'+'\t'.join(Z)+'\n')

open(OUTDIR+'/header.txt','w').write('\t'.join(header)+'\n')

jobs=[]
for edge in data.PCC_POOL:
    ps=edge.split(':')
    p1=ps[0]
    p2=ps[1]
    try:
        p1_old_exp=[]#data.EXP[p1]
        p2_old_exp=[]#data.EXP[p2]
        i=0
        while i<data.POOL_LENGTH:
            if data.EXP[p1][i] * data.EXP[p2][i] !=0:
                p1_old_exp.append(data.EXP[p1][i])
                p2_old_exp.append(data.EXP[p2][i])
            i+=1
        p1_this_list=data.EXP[p1]
        p2_this_list=data.EXP[p2]
        pcc_old = data.PCC_POOL[edge][0]
        pcc_length = data.PCC_POOL[edge][1]
        p= multiprocessing.Process(target=SINGLE, args=(p1,p2,p1_old_exp,p2_old_exp,p1_this_list, p2_this_list, pcc_old, pcc_length))
        p.start()
        jobs.append(p)
    except Exception as e:
        print e
    if len(jobs)>=PROC_LIMIT:
        for p in jobs:
            p.join()
        jobs=[]
for p in jobs:
    p.join()

#subprocess.Popen('cat '+OUTDIR+'/*.ssn > '+OUTDIR+'/DATA.SSN',shell=True).wait()
#subprocess.Popen('cat '+OUTDIR+'/header.txt '+OUTDIR+'/DATA.SSN > '+OUT_FILE,shell=True).wait()
#################################    
FILES=os.listdir(OUTDIR)
fo=open(OUT_FILE,'w')
fo.write(open(OUTDIR+'/header.txt').read())
for one in FILES:
    if '.ssn' in one: 
        fo.write(open(OUTDIR+'/'+one).read())
################################# 




