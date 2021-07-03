from __future__ import print_function
from builtins import range
import numpy as np
import os
from astropy_healpix import HEALPix
from astropy import units as u
import sys
from astropy.table import Table

hp = HEALPix(nside=256)
print(hp.npix,'total healpix pixels on sky')
area=hp.pixel_area.value*3283
print('area of one healpix is',area,'sq. deg')
print('typical source separation should be <',np.sqrt(area)*3600.0,'arcsec')

for f in sys.argv[1:]:
    path=f.split('/')[0]
    t=Table.read(f)
    if 'ID_ra' in t.colnames:
        ra='ID_ra'
        dec='ID_dec'
    else:
        ra='ID_RA'
        dec='ID_DEC'
    nsources=len(t)
    nid=np.sum(~np.isnan(t[ra]))
    nz=np.sum(~np.isnan(t['z_best']))
    
    pixels=np.unique(hp.lonlat_to_healpix(t[ra].quantity.value*u.deg,t[dec].quantity.value*u.deg))

    #print(len(pixels),'total distinct pixels')
    #print('Area is',len(pixels)*area,'sq. deg')
    print("%-15s %6i %6i %4.1f %6i %4.1f %6.1f" % (path,nsources,nid,nid*100.0/nsources,nz,nz*100.0/nsources,len(pixels)*area))
