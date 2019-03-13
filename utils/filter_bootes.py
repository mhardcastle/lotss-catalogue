#!/usr/bin/python

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import fits
import sys
from astropy.wcs import WCS
from read_maps import ReadMaps

# arguments: input cat, output cat

t=Table.read(sys.argv[1])
spitzer=fits.open('/beegfs/lofar/deepfields/Bootes_optical/SDWFS/I2_bootes.v32.fits')
iband=ReadMaps('/beegfs/lofar/deepfields/Bootes_optical/NDWFS/I/*_con.fits')

w=WCS(spitzer[0].header)

pos=w.wcs_world2pix(t['RA'],t['DEC'],0)
print len(pos)

filter=[]
for i,r in enumerate(t):

    x=int(pos[0][i])
    y=int(pos[1][i])
    try:
        inmask=(spitzer[0].data[y,x]!=0)
    except:
        inmask=False
    if inmask:
        # now look up position in NDWFS images
        res=iband.find_pos(r['RA'],r['DEC'])
        if res is None:
            inmask=False
    filter.append(inmask)
    print i,x,y,inmask

t_new=t[filter]
t_new.write(sys.argv[2],overwrite=True)
