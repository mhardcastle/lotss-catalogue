# merge hp maps

from astropy.table import Table
import numpy as np
from tqdm import tqdm
import healpy as hp
from astropy.io import fits
import os

NSIDE=8192 # was 4096

import json

with open('/data/lofar/DR3/vlow_sub_healpix_map.json') as infile:
    j=json.load(infile)

pixels=hp.pixelfunc.nside2npix(NSIDE)

temperature=np.ones(pixels)*np.nan

beamarea=18.12944056730881*(0.004166666666666667*np.pi/180)**2 # sr

for k in tqdm(j):
    f=f'/data/lofar/DR3/healpix_mosaics/{k}/vlow-hptable-{NSIDE}.fits'
    if not os.path.isfile(f): continue
    t=Table.read(f)
    # Unit conversions
    # The median flux here is in Jy/beam
    # I_nu is flux/area of beam in sr
    # T = I_nu c^2/2k nu^2
    I_nu=t['Flux per pixel']*1e-26/beamarea
    temps=I_nu * 3e8**2/(2*1.38e-23*144e6**2)
    temperature[t['PIXEL']]=temps
    if np.max(np.abs(temps))>1e5:
        print(k,np.max(temps))
    
hp.fitsfunc.write_map(f'/data/lofar/mjh/temperature_hp-{NSIDE}.fits',temperature,overwrite=True)
