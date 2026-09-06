# convert mosaic to healpix
# started life as mos2hpmap_hpm.py but this version will make a temperature map for the VLOW feathering

import sys
import os
import healpy as hp
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.table import Table
import numpy as np
from tqdm import tqdm
from scipy import ndimage
import json
from astropy import units as u
from astropy.coordinates import SkyCoord

def subpixels_in_pixel(ipix, nside_1, nside_2):
    """
    ChatGPT
    Return all RING-ordered HEALPix pixels at nside_2 contained inside
    the RING-ordered pixel ipix at nside_1.

    Parameters
    ----------
    ipix : int
        Pixel number at nside_1 (RING ordering).
    nside_1 : int
        Coarse NSIDE.
    nside_2 : int
        Fine NSIDE (must be > nside_1).

    Returns
    -------
    numpy.ndarray
        Array of pixel numbers at nside_2 (RING ordering).
    """

    if nside_2 <= nside_1:
        raise ValueError("nside_2 must be larger than nside_1")

    # Convert input pixel from RING -> NEST
    ipix_nest = hp.ring2nest(nside_1, ipix)

    # Ratio between resolutions
    factor = nside_2 // nside_1

    if nside_1 * factor != nside_2:
        raise ValueError("nside_2 must be an integer multiple of nside_1")

    # Number of subpixels per parent pixel
    n_subpix = factor**2

    # In NEST ordering, children are contiguous
    first_child = ipix_nest * n_subpix
    children_nest = np.arange(first_child, first_child + n_subpix)

    # Convert back to RING ordering
    children_ring = hp.nest2ring(nside_2, children_nest)

    return sorted(children_ring)

if __name__=='__main__':

    NSIDE=8192 # was 4096
    NSIDE_BIG=16
    with open('/data/lofar/DR3/vlow_sub_healpix_map.json') as infile:
        j=json.load(infile)
    for bp in tqdm(j):
        big_pixel=int(bp)
        #pixels=subpixels_in_pixel(big_pixel,NSIDE_BIG,NSIDE)
        wd=f'/data/lofar/DR3/healpix_mosaics/{big_pixel}'
        if not os.path.isdir(wd): continue
        os.chdir(wd)
        if os.path.isfile(f'vlow-hptable-{NSIDE}.fits'): continue
        hdu=fits.open('vlow-sub-mosaic-blanked.fits')
        #print(hdu[0].data.shape)
        #print(np.sum(~np.isnan(hdu[0].data)),'non-blanked pixels')
        wcs=WCS(hdu[0].header)
        y=np.arange(0,hdu[0].data.shape[0])
        x=np.arange(0,hdu[0].data.shape[1])
        xx,yy=np.meshgrid(x,y)
        ra,dec=wcs.all_pix2world(xx,yy,0)
        sc=SkyCoord(ra*u.deg,dec*u.deg,frame='icrs')
        gal=sc.galactic
        l=gal.l.value
        b=gal.b.value
        hpm=hp.pixelfunc.ang2pix(NSIDE,l,b,lonlat=True)
        pixels=np.unique(hpm[~np.isnan(hdu[0].data)])
        """
        counts=np.zeros_like(pixels)
        fluxes=np.ones_like(pixels)*np.nan
        for i,p in enumerate(tqdm(pixels)):
            mask=(hpm==p)
            masked=hdu[0].data[mask]
            counts[i]=np.sum(~np.isnan(masked))
            fluxes[i]=np.nanmedian(masked)
        """
        labels = hpm.ravel()
        values = hdu[0].data.ravel()

        # Keep only finite values
        valid = ~np.isnan(values)
        labels_valid = labels[valid]
        values_valid = values[valid]
        if len(values_valid)==0: continue

        # ---- Counts ----
        #counts_dict = np.bincount(labels_valid)

        # Extract counts for requested pixels
        #counts = counts_dict[pixels]

        # ---- Medians ----
        # ndimage.median works by label
        fluxes = ndimage.median(
            values_valid,
            labels=labels_valid,
            index=pixels
        )

        t=Table([pixels,fluxes],names=['PIXEL','Flux per pixel'])
        t.write(f'vlow-hptable-{NSIDE}.fits',overwrite=True)



