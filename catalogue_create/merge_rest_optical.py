# Merge Ken's restframe optical catalogue
import numpy as np
from astropy.table import Table,hstack

t=Table.read('LOFAR_HBA_T1_DR1_merge_ID_optical_v1.0.fits')
t2=Table.read('../../LOFAR_HBA_T1_DR1_merge_ID_optical_v1.0_restframe_mags_v0.3.fits')
del(t2['Source_Name'])
del(t2['z_best'])
tm=hstack((t,t2))
restcols=[n for n in tm.colnames if 'rest' in n]
for c in restcols:
    tm[c]=np.where(tm[c]==-99,np.nan,tm[c])

tm.write('LOFAR_HBA_T1_DR1_merge_ID_optical_v1.0_restframe.fits',overwrite=True)
