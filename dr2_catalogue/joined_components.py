from __future__ import print_function
from astropy.table import Table,vstack
import glob

fields=['Spring-40-45','Spring-60-65','Winter','Fall']
versions=['v0.4','v0.8','v0.5','v0.6'] # matches version 0.1 source cat
tables=[]

for f,v in zip(fields,versions):
    g=f+'/components-'+v+'.fits'
    print(g)
    tables.append(Table.read(g))
    tables[-1]['Field']=f

t=vstack(tables)

t.write('combined_components_v0.1.fits',overwrite=True)

