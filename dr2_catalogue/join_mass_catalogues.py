# Create one mass catalogue from the two we have

from __future__ import print_function
import glob
from astropy.table import Table,hstack,MaskedColumn,join
import numpy as np
from astropy.io import fits

kdcat1='/beegfs/lofar/duncan/combined_release_v1.0_opt_mass_cols.fits'
kdcat2='/beegfs/lofar/duncan/combined_release_v1.0b_zspec_mass_cols.fits'

tkd1=Table.read(kdcat1)
tkd2=Table.read(kdcat2)
t=join(tkd1,tkd2,keys='Source_Name',join_type='left')

selection=~t['K_rest_2'].mask
for c in list(t.colnames):
    if c.endswith('_1'):
        col=c.replace('_1','')
        t[col]=np.where(selection,t[col+'_2'],t[col+'_1'])

t['Mass_from']=np.where(selection,'Spec','Phot')

for c in list(t.colnames):
    if c.endswith('_1') or c.endswith('_2'):
        del(t[c])

t['Mass_from']=np.where(np.isnan(t['Mass_median']),'',t['Mass_from'])
        
t.write('fixed_mass_table.fits',overwrite=True)
