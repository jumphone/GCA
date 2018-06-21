fi=open('STEM.txt')
stem=set()
for line in fi:
    seq=line.rstrip().split('\t')
    stem.add(seq[0])

fi.close()
fi=open('trrust_rawdata.human.tsv')
fo=open('trrust_rawdata.human.tsv.tgstem','w')
for line in fi:
    seq=line.rstrip().split('\t')
    if seq[1] in stem:
        fo.write(line)

