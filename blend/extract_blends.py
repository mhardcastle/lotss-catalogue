from astropy.table import Table
import numpy as np

t=Table.read('LOFAR_HBA_T1_DR1_merge_ID_v1.0.fits')

t=t[((t['ID_flag']>=61) & (t['ID_flag']<=63))]

t.write('blends.fits',overwrite=True)
