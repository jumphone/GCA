fi=open('MODE.summary.txt')
fo=open('MODE_MAT.txt','w')


TFTG={}

for line in fi:
    input_file=line.rstrip()
    tag='TF_'+input_file.split('TF_')[1].split('.txt')[0]
    tag=tag.replace('_','.')
    #print tag
    ftmp=open(input_file)
    ftmp.readline()
    out=[tag]
    header=[]
    for l in ftmp:
        sss=l.rstrip().split('\t') 
        header.append(sss[0])
        if sss[-1]=='TG_HI;TF_POS;':
            out.append('1')
        else:
            out.append('0') 
    TFTG[tag]='\t'.join(out)
 

fo.write('\t'.join(header)+'\n')
for tag in TFTG:
    fo.write(TFTG[tag]+'\n')
