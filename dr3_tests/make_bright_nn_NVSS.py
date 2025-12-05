# Generate a bright cut catalogue with a nearest neighbour filter

from astropy.table import Table
import healpy as hp
import numpy as np
from tqdm import tqdm

def separation(c_ra,c_dec,ra,dec):
    # all values in degrees
    return np.sqrt((np.cos(c_dec*np.pi/180.0)*(ra-c_ra))**2.0+(dec-c_dec)**2.0)

def filter_catalogue(t,c_ra,c_dec,radius):
#    r=np.sqrt((np.cos(c_dec*np.pi/180.0)*(t['RA']-c_ra))**2.0+(t['DEC']-c_dec)**2.0)
    r=separation(c_ra,c_dec,t['RA'],t['DEC'])
    return t[r<radius]

def select_isolated_sources(t,radius,progress=False):
    t['NN_dist']=np.nan
    if progress:
        iterable=tqdm(t)
    else:
        iterable=t
    for r in iterable:
        filt=t['pixel']==r['pixel']
        dist=3600.0*separation(r['RA'],r['DEC'],t['RA'][filt],t['DEC'][filt])
        # dist=np.sqrt((np.cos(c_dec*np.pi/180.0)*(t['RA']-r['RA']))**2.0+(t['DEC']-r['DEC'])**2.0)*3600.0
        dist.sort()
        if len(dist)>1:
            r['NN_dist']=dist[1]

    t=t[t['NN_dist']>radius]
    return t


fluxcut=0.01 # Jy
sep=120 # arcsec
NSIDE=64

t=Table.read('/beegfs/general/mjh/bootstrap/NVSS.fits')

t=t[t['Total_flux']>fluxcut]
print(len(t))

pixels=hp.ang2pix(NSIDE,t['RA'],t['DEC'],lonlat=True)
t['pixel']=pixels

t=select_isolated_sources(t,sep,progress=True)

t.write('NVSS-bright-cut-NN1.fits')

