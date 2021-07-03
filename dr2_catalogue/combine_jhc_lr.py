from __future__ import print_function
from astropy.table import Table
import numpy as np

t=Table.read('/beegfs/lofar/mjh/rgz/Spring-60-65/sources-v0.3.fits')
t2=Table.read('/beegfs/lofar/jcroston/surveys/dr2_hosts/ww/bbww_tests/13h_RLC_10mJy_15arcsec_LR1.fits')

skip=0
replace=0
notfound=0
for r in t2:
    search=t['Source_Name']==r['Source_Name']
    if not np.any(search):
        print(r['Source_Name'],'not found!')
        notfound+=1
        continue
    index=np.argmax(search)
    if t[index]['Position_from']=='Visual inspection' and not np.isnan(t[index]['optRA']):
        print(r['Source_Name'],'skipped as visual inspection position exists',index)
        skip+=1
    else:
        print(r['Source_Name'],'replacing position derived from',t[index]['Position_from'],index)
        replace+=1

print(skip,replace,notfound,skip+replace)
        
