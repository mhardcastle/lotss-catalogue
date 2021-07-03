from __future__ import print_function
import glob
import os
from astropy.table import Table
import sys
import numpy as np

dir=os.getcwd()
field=os.path.basename(dir)
print('field is',field)
optcat='optical.fits'

g=sorted(glob.glob('sources-v*.fits'))

infile=g[-1]
print('Processing',infile)

mergeout=infile.replace('sources','merged')

os.system('/soft/topcat/stilts tskymatch2 ra1=optRA dec1=optDec ra2=RA dec2=DEC error=2.0 join=all1 in1=%s in2=%s out=%s' % (infile,optcat,mergeout))

t=Table.read(mergeout)
finalname=mergeout.replace('merged','final')
print('Remove unnecessary columns')

cols=['UID_L_1','NoID','Art_prob','Hostbroken_prob','Imagemissing_prob','Zoom_prob']

for c in cols:
    if c in t.colnames:
        print(c,end=' ')
        sys.stdout.flush()
        t.remove_columns(c)

print('\nRename columns:',end=' ')

for oldcol,column in [('RA_1','RA'),('DEC_1','DEC'),('RA_2','ID_RA'),('DEC_2','ID_DEC'),('ID','ID_NAME'),('UID_L_2','UID_L')]:

    if oldcol in t.colnames:
        t[oldcol].name=column
        print(oldcol,'->',column,end=': ')
        sys.stdout.flush()

print('\nRemove whitespace padding:',end='')
sys.stdout.flush()
stringcols=[n for (n,ty) in t.dtype.descr if 'S' in ty]
for c in stringcols:
    print(c,end=' ')
    sys.stdout.flush()
    t[c]=[s.rstrip() for s in t[c]]
print()

print('Sorting')
t.sort('RA')

print('Writing %s to disk' % finalname)
t.write(finalname,overwrite=True)
