fi=open('MODE_MAT.txt')
fo=open('MODE_MAT.txt.tag','w')

header=fi.readline().rstrip().split('\t')
for one in header:
    tag=one.rstrip().split('_')[0].replace('X','')
    if 'MGH' not in tag:
        tag='MGH'+tag
    fo.write(tag+'\n')




