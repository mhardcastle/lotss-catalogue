from astropy.table import Table
import numpy as np

t = Table.read('LOFAR_HBA_T1_DR1_merge_ID_optical_f_v1.2_restframe.fits')

mlist = np.unique(t['Mosaic_ID'])


for tt in mlist:
    matched = []
    if 'Hetde' in tt:
        for m in mlist:
            if tt != m:
                if tt in m:
                    matched.append(m)
    print tt, matched
    if len(matched) == 1 :
        t['Mosaic_ID'][t['Mosaic_ID']==tt] = matched[0]
    elif len(matched) > 1:
        print 'error'
