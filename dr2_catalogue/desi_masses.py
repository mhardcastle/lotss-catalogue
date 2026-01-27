# dirty merge of the DESI masses where we correct if the redshift is close enough

import numpy as np
import glob
import sys
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.table import MaskedColumn, Table, join, vstack, hstack

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

t1=Table.read('combined-release-v1.2-LM.fits')
t2=Table.read('combined-release-v1.1-LM_opt_mass.fits')

z_best_old=t2['z_best']
zchange=t1['z_best']!=z_best_old 

print(np.sum(zchange),'redshifts have changed')

# figure of merit -- small values new z is very close to old
fom=np.abs(z_best_old-t1['z_best'])/(1+t1['z_best'])

ldr=cosmo.luminosity_distance(t1['z_best'])/cosmo.luminosity_distance(z_best_old)

correct=zchange & (fom<0.05)
blank=zchange & (fom>=0.05)

print(np.sum(correct),'masses are corrected')
print(np.sum(blank),'masses are blanked')
      
masscols=[k for k in t2.colnames if k.startswith('Mass_')]
dsm=u.dex(u.M_sun)
for k in masscols:
    print(k)
    t1[k]=np.array(t2[k])
    t1[k][correct]+=2*np.log10(ldr[correct])
    t1[k][blank]=np.nan
    t1[k]*=dsm
    

#t['badmass']=z_best_old!=t['z_best']

t1.write('combined-release-v1.2-LM_scaledmass.fits',overwrite=True)
