# convert mosaic to healpix

import sys
import os
import healpy as hp
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.table import Table
import numpy as np
from tqdm import tqdm

NSIDE=4096
NSIDE_BIG=16

big_pixel=
os.chdir('/beegfs/lofar/DR3/healpix_mosaics/'+sys.argv[1])
hdu=fits.open('mosaic-blanked.fits')
rms=fits.open('mosaic-blanked--final.rms.fits')
print(hdu[0].data.shape)
rmsdata=rms[0].data[0,0]
wcs=WCS(hdu[0].header)
y=np.arange(0,hdu[0].data.shape[0])
x=np.arange(0,hdu[0].data.shape[1])
xx,yy=np.meshgrid(x,y)
ra,dec=wcs.wcs_pix2world(yy,xx,0)
hpm=hp.pixelfunc.ang2pix(NSIDE,ra,dec,lonlat=True)
pixels=sorted(np.unique(hpm))
counts=np.zeros_like(pixels)
fluxes=np.ones_like(pixels)*np.nan
rmsv=np.ones_like(pixels)*np.nan
for i,p in enumerate(tqdm(pixels)):
    mask=(hpm==p)
    counts[i]=np.sum(~np.isnan(hdu[0].data[mask]))
    if counts[i]>0:
        fluxes[i]=np.nansum(hdu[0].data[mask])/counts[i]
        rmsv[i]=np.nanmean(rmsdata[mask])

t=Table([pixels,counts,fluxes,rmsv],names=['PIXEL','Count','Flux per pixel','RMS'])
t.write('hptable.fits',overwrite=True)

