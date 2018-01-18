# Search the zoom over-ride files for sources marked as artefacts and make a list for use in Wendy's code

import glob

outfile=open('artefacts.txt','w')

dir='/data/lofar/mjh/hetdex_v4/zoom/'
g=glob.glob(dir+'ILT*.txt')
for f in g:
    if 'table-list' in f:
        continue
    n=f.split('/')[-1].replace('.txt','')

    lines=open(f).readlines()
    if lines[0].rstrip()=='## Deleted':
        if len(lines)==1:
            outfile.write(n+'\n')
        else:
            outfile.write(''.join(lines[1:]))
            
