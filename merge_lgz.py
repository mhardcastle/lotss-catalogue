## This code is NOT used for the merge because it gives slightly strange results.
## I don't know why but running topcat/stilts is easy...

from astropy.table import Table
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

ocat=Table.read('/data/lofar/mjh/hetdex_ps1_allwise_photoz_v0.2.fits')
t=Table.read('LOFAR_HBA_T1_DR1_merge_ID_v0.5.fits')
opcat=SkyCoord(ocat['ra']*u.deg, ocat['dec']*u.deg)

indices=np.where(~np.isnan(t['ID_ra']))[0]

withid=t[indices]
withidc=SkyCoord(withid['RA'], withid['DEC'])
idx, d2d, d3d = withidc.match_to_catalog_sky(opcat)

# now we have a set of indices into ocat
# to do a left join we make a dummy column to join on

ocat['Index'][idx]=indices+1
t['Index']=0
t['Index']=range(1,1+len(t))

j=join(t,ocat,keys='Index',join_type='left')
