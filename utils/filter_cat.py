#!/usr/bin/python

from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io import fits
import sys
from astropy.wcs import WCS

# arguments: input cat, mask fits, output cat

t=Table.read(sys.argv[1])
mask=fits.open(sys.argv[2])
w=WCS(mask[0].header)

pos=w.wcs_world2pix(t['RA'],t['DEC'],0)
print len(pos)

filter=[]
for i,r in enumerate(t):

    x=int(pos[0][i])
    y=int(pos[1][i])
    try:
        inmask=(mask[0].data[y,x]>0)
    except:
        inmask=False
    filter.append(inmask)
    print i,x,y,inmask

t_new=t[filter]
t_new.write(sys.argv[3],overwrite=True)
