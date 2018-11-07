#!/usr/bin/python

from astropy.table import Table
import numpy as np

t=Table.read('../blend/merge_out.fits')

print 'Stripping...'
t['Source_Name']=[s.rstrip() for s in t['Source_Name']]
    
tr=Table.read('tmissing_psmatch.fits')

print 'Merging'
for r in tr:
    id=r['Source_Name'].rstrip()
    mask=t['Source_Name']==id
    i=np.argmax(mask)
    print i,id
    t[i]['ID_flag']=22
    t[i]['ID_ra']=r['ra_2']
    t[i]['ID_dec']=r['dec_2']

print 'Fixing up 2MASS names'
for i,r in enumerate(t):
    if '2MASS' in r['ID_name']:
        t[i]['ID_flag']=22
        t[i]['ID_name']=r['ID_name'].replace('2MASS ','2MASX J')
        print i, r['ID_name'],t[i]['ID_name']
        
print 'Writing to disk'
t.write('merge_out_fixed.fits',overwrite=True)
