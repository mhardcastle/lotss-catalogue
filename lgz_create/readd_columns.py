from astropy.table import Table
import numpy as np

t=Table.read('LGZ_philip.fits')
t_full=Table.read('/data/lofar/wwilliams/hetdex/LoTSS-DR1-July21-2017/LOFAR_HBA_T1_DR1_merge_ID_v0.7.fits')
#for i in range(len(t_full)):
#    t_full[i]['Source_Name']=t_full[i]['Source_Name'].rstrip()

cols=['Maj','Min','PA']
for c in cols:
    t.remove_column(c)
    t[c]=np.nan

for i,r in enumerate(t):
    n=r['Source_Name']
    print i,n
    try:
        r2=t_full[t_full['Source_Name']==n][0]
        if (r2['ID_flag']==5):
            for c in cols:
                t[i][c]=r2[c]
    except IndexError:
        print 'Source not found'
    
t=t[~np.isnan(t['Maj'])]
t.write('LGZ_philip_fixed.fits',overwrite=True)
