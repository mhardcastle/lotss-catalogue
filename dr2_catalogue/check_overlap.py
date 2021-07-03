from __future__ import print_function
import numpy as np
import os
from astropy_healpix import HEALPix
from astropy import units as u
import sys
from astropy.table import Table,vstack

t1=Table.read('Winter/final-v0.5_photoz_v0.1_joined-physical.fits')
t2=Table.read('Spring-40-45/final-v0.4_photoz_v0.1_north-physical.fits')
t3=Table.read('Spring-60-65/final-v0.8_photoz_v0.1-physical.fits')
t4=Table.read('Fall/final-v0.6_photoz_v0.1-physical.fits')
t1['field']='8h'
t2['field']='13h40'
t3['field']='13h60'
t4['field']='0h'

hp = HEALPix(nside=256)

t1_pixel=np.array(hp.lonlat_to_healpix(t1['RA'].quantity.value*u.deg,t1['DEC'].quantity.value*u.deg))
t2_pixel=np.array(hp.lonlat_to_healpix(t2['RA'].quantity.value*u.deg,t2['DEC'].quantity.value*u.deg))

t1_pixels=np.unique(t1_pixel)
t2_pixels=np.unique(t2_pixel)

overlap=np.intersect1d(t1_pixels,t2_pixels)

from_t1=np.setdiff1d(t1_pixels,overlap) # the ones not in overlap
from_t2=np.setdiff1d(t2_pixels,overlap) 

print(len(overlap),'overlapping pixels')
for p in overlap:
    if np.sum(t1_pixel==p)>np.sum(t2_pixel==p):
        print('Selecting t1 for pixel',p)
        from_t1=np.append(from_t1,[p])
    else:
        print('Selecting t2 for pixel',p)
        from_t2=np.append(from_t2,[p])

t1filter=np.zeros_like(t1,dtype=np.bool)
for p in from_t1:
    t1filter|=(t1_pixel==p)

t2filter=np.zeros_like(t2,dtype=np.bool)
for p in from_t2:
    t2filter|=(t2_pixel==p)
    
combined=vstack([t1[t1filter],t2[t2filter],t3,t4])
combined.sort('RA')
combined.write('combined-release-v0.1.fits',overwrite=True)
