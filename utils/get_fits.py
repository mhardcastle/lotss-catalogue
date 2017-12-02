#!/usr/bin/python

import sys
import os
from subim import extract_subim
from download_image_files import LofarMaps
from astropy.coordinates import SkyCoord
import astropy.units as u

if len(sys.argv)==3:
    ra=float(sys.argv[1])
    dec=float(sys.argv[2])
    sc=SkyCoord(ra*u.deg,dec*u.deg,frame='icrs')
    s=sc.to_string(style='hmsdms',sep='',precision=2)
    name=(str('ILTJ'+s).replace(' ','')[:-1])
elif len(sys.argv)==2:
    s=sys.argv[1][4:]
    coord=s[0:2]+':'+s[2:4]+':'+s[4:9]+' '+s[9:12]+':'+s[12:14]+':'+s[14:]
    sc = SkyCoord(coord,unit=(u.hourangle,u.deg))
    ra=sc.ra.value
    dec=sc.dec.value
    print 'Parsed coordinates to ra=%f, dec=%f' % (ra,dec)
    name=sys.argv[1]
else:
    print 'Call me with the name of a source or RA, Dec in degrees'
    

lm=LofarMaps()
filename=os.environ['IMAGEDIR']+'/'+lm.find(ra,dec)
hdu=extract_subim(filename,ra,dec,0.25)
hdu.writeto(name+'.fits',clobber=True)
