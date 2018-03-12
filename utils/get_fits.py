#!/usr/bin/python

import sys
import os
from subim import extract_subim
from download_image_files import LofarMaps
from astropy.coordinates import SkyCoord
import astropy.units as u

def save_fits(ra,dec,name=None,size=0.25):
    if name is None:
        sc=SkyCoord(ra*u.deg,dec*u.deg,frame='icrs')
        s=sc.to_string(style='hmsdms',sep='',precision=2)
        name=(str('ILTJ'+s).replace(' ','')[:-1])
    lm=LofarMaps()
    filename=os.environ['IMAGEDIR']+'/'+lm.find(ra,dec)
    hdu=extract_subim(filename,ra,dec,size)
    outname=name+'.fits'
    hdu.writeto(outname,clobber=True)
    return outname

if __name__=='__main__':
    
    if len(sys.argv)==3:
        ra=float(sys.argv[1])
        dec=float(sys.argv[2])
        save_fits(ra,dec)
    elif len(sys.argv)==4:
        ra=float(sys.argv[1])
        dec=float(sys.argv[2])
        size=float(sys.argv[3])
        save_fits(ra,dec,size=size)
    elif len(sys.argv)==2:
        s=sys.argv[1][4:]
        coord=s[0:2]+':'+s[2:4]+':'+s[4:9]+' '+s[9:12]+':'+s[12:14]+':'+s[14:]
        sc = SkyCoord(coord,unit=(u.hourangle,u.deg))
        ra=sc.ra.value
        dec=sc.dec.value
        print 'Parsed coordinates to ra=%f, dec=%f' % (ra,dec)
        name=sys.argv[1]
        save_fits(ra,dec,name=name)
    else:
        print 'Call me with the name of a source or RA, Dec in degrees'
