# Merge Ken's restframe optical catalogue
import numpy as np
from astropy.table import Table,hstack

t=Table.read('LOFAR_HBA_T1_DR1_merge_ID_optical_f_v1.2.fits')
t2=Table.read('LOFAR_HBA_T1_DR1_merge_ID_optical_f_v1.2_restframe_mags.fits')
del(t2['Source_Name'])
del(t2['z_best'])
tm=hstack((t,t2))
restcols=[n for n in tm.colnames if 'rest' in n]
for c in restcols:
    tm[c]=np.where(tm[c]==-99,np.nan,tm[c])

tm.write('LOFAR_HBA_T1_DR1_merge_ID_optical_f_v1.2_restframe.fits',overwrite=True)
