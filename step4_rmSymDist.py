################################################
#Author: Feng Zhang
#Email: 15110700005@fudan.edu.cn
#Date: 2018.06.12
################################################

import sys
from scipy import stats
import statsmodels.stats.multitest as stm

import numpy as np

print'''
$1 INPUT
$2 OUTPUT
'''


INPUT = sys.argv[1]
OUTPUT = sys.argv[2]

outlierConstant=1.5
#def removeOutliers(x, outlierConstant):
def removeOutliers(x):
    a = np.array(x)
    upper_quartile = np.percentile(a, 75)
    lower_quartile = np.percentile(a, 25)
    IQR = (upper_quartile - lower_quartile) * outlierConstant
    quartileSet = (lower_quartile - IQR, upper_quartile + IQR)
    resultList = []
    for y in a.tolist():
        if y >= quartileSet[0] and y <= quartileSet[1]:
            resultList.append(y)
    return resultList



fi=open(INPUT)
fo=open(OUTPUT+'.tmp','w')
fi.readline()
OUT=[]
P=[]
for line in fi:
    seq=line.rstrip().split('\t')
    tag=seq[0]
    Z_LIST=[]
    NEG=[]
    POS=[]
    for one in seq[1:]:
        if one !='NA':
            Z_LIST.append(float(one))
    Z_MED = np.median(Z_LIST)
    Z_LIST=removeOutliers(Z_LIST)
    Z_OUT=[]
    for one in Z_LIST:
            Z_OUT.append(str(one))
            value=float(one)-Z_MED
            if value<0:
                NEG.append(abs(value))
            elif value>0:
                POS.append(abs(value))
    st, pv=stats.ks_2samp(NEG,POS)
    #OUT.append([pv,tag])
    P.append(pv)
    OUT.append([tag,str(pv),','.join(Z_OUT)])


adP=stm.multipletests(P,method='fdr_bh')[1]
#adP=stm.multipletests(P,method='bonferroni')[1]
i=0
SIG=set()
while i<len(OUT):
    flag='_NO_'
    if adP[i]<0.05:
        flag='_SIG_'
        SIG.add(OUT[i][0])
    #fo.write(OUT[i][0]+'\t'+OUT[i][1]+'\t'+str(adP[i])+'\t'+flag+'\t'+OUT[i][2]+'\n')
    fo.write(OUT[i][0]+'\t'+OUT[i][1]+'\t'+str(adP[i])+'\t'+flag+'\n')
    i+=1

fi.close()
fi=open(INPUT)
fo=open(OUTPUT,'w')
fo.write(fi.readline())
for line in fi:
    seq=line.rstrip().split('\t')
    if seq[0] in SIG:
        fo.write(line)
fo.close()
fi.close()

