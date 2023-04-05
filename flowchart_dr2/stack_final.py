from astropy.table import Table,Column, vstack
import numpy as np

t0 = Table.read('LoTSS_DR2_v110.srl_0h.lr-full.sorted_step3_flux4.fits')
t13 = Table.read('LoTSS_DR2_v110.srl_13h.lr-full.sorted_step3_flux4.fits')

# add missing columns from t13 to t0
# these are all boolean zeros
new0 = ['WEAVE_priority1a','WEAVE_priority2','WEAVE_priority3','HETDEX']
for n in new0:
    t0.add_column(Column(data=np.zeros(len(t0),dtype=bool),name=n))

t = vstack([t0,t13])

t.write('LoTSS_DR2_v110.srl.lr-full.sorted_step3_flux4.fits', overwrite=True)
