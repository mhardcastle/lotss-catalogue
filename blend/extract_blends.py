from astropy.table import Table
import numpy as np

t=Table.read('LOFAR_HBA_T1_DR1_merge_ID_v1.2.fits')

t=t[((t['ID_flag']>=41) & (t['ID_flag']<=42))]

t.write('blends.fits',overwrite=True)
