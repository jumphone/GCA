fi=open('MODE.summary.txt')
fo=open('MODE_MAT.txt','w')


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


