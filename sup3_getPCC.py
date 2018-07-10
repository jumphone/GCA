import cPickle as pickle
import sys
from scipy import stats

INPUT=sys.argv[1]
OUTPUT=sys.argv[2]


class Person_PCC_Data:
    def __init__(self,  EXP , POOL_LENGTH , PCC_POOL ):
        self.EXP = EXP
        self.POOL_LENGTH = POOL_LENGTH
        self.PCC_POOL = PCC_POOL

INDEX_DATA_FILE=INPUT

print 'loading...'
fdata = open(INDEX_DATA_FILE)
data = pickle.load(fdata)
fdata.close()
print "loading done !"
fo=open(OUTPUT,'w')
for pair in data.PCC_POOL:
    fo.write(pair.replace(':','.And.')+'\t'+str(data.PCC_POOL[pair][0])+'\t'+str(data.PCC_POOL[pair][1])+'\n')
