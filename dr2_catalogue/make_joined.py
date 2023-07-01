from __future__ import print_function
import numpy as np
import os
import sys
from astropy.table import Table,vstack
from surveys_db import SurveysDB

# simplified from check_overlap.py since no duplicated sky area
# merge in legacy coverage and component tables

fields=['Spring','Fall']
versions=['v2.3','v2.5']
pzname=['v0.3_joined','v0.3_south']
outversion='v1.1'
tables=[]

print('reading tables')
for f,v,p in zip(fields,versions,pzname):
    g=f+'/final-'+v+'_photoz_'+p+'-physical.fits'
    print(g)
    tables.append(Table.read(g))
    tables[-1]['Field']=f

print('stacking and sorting')
t=vstack(tables)
t.sort('RA')

print('temp fix for ID_flag')
t['ID_flag']=np.where(t['Position_from']=='Ridge line code',13,t['ID_flag'])

print('add legacy coverage flag')

with SurveysDB() as sdb:
    sdb.execute('select id,gz_status from fields where dr2')
    res=sdb.cur.fetchall()

s={}
for r in res:
    s[r['id']]=r['gz_status']
    
t['Legacy_Coverage']=np.zeros_like(t,dtype=bool)
flag=[]

for k in s:
    if s[k]!='Not usable':
        t['Legacy_Coverage']|=(t['Mosaic_ID']==k)

print('write')
t.write('combined-release-'+outversion+'.fits',overwrite=True)

print('read in component tables')
tables=[]

for f,v in zip(fields,versions):
    g=f+'/components-'+v+'.fits'
    print(g)
    tables.append(Table.read(g))
    tables[-1]['Field']=f

print('stack and sort')
t=vstack(tables)
t.sort('RA')

print('write')
t.write('combined-components-'+outversion+'.fits',overwrite=True)
