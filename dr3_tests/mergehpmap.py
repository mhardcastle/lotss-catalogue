# merge hp maps

import glob
from astropy.table import Table
import numpy as np
from tqdm import tqdm
import healpy as hp
from astropy.io import fits

NSIDE=4096

g=glob.glob('/beegfs/lofar/DR3/mosaics-new/*/hptable.fits')

pixels=hp.pixelfunc.nside2npix(NSIDE)

flux=np.ones(pixels)*np.nan
rms=np.ones(pixels)*np.nan
area=np.ones(pixels)*np.nan

for f in tqdm(g):
    t=Table.read(f)
    hdu=fits.open(f.replace('hptable.fits','mosaic-blanked.fits'))
    if hdu[0].header['BMAJ']>0.0023:
        t['Flux per pixel']/=1.5**2
    filt=t['Count']>0
    flux[t['PIXEL'][filt]]=t['Flux per pixel'][filt]
    rms[t['PIXEL'][filt]]=t['RMS'][filt]
    area[t['PIXEL'][filt]]=np.where(np.isnan(area[t['PIXEL'][filt]]),t['Count'][filt],area[t['PIXEL'][filt]]+t['Count'][filt])

t=Table([flux,rms,area],names=["Flux","RMS","area"])
t.write(f'fulltable-{NSIDE}.fits',overwrite=True)
    
hp.fitsfunc.write_map('flux_hp.fits',flux,overwrite=True)
hp.fitsfunc.write_map('area_hp.fits',area,overwrite=True)
hp.fitsfunc.write_map('rms_hp.fits',rms,overwrite=True)
