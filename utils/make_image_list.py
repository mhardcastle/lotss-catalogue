# Like download_image_files, but assumes that the relevant files are
# already there, and so just finds the nearest map from a list.

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
import os
import sys
import glob
import numpy as np
import warnings
from astropy.utils.exceptions import AstropyWarning

class MapsList(object):
    def __init__(self,globstring,hdu_no=0):
        warnings.simplefilter('ignore',category=AstropyWarning)
        print 'Ingesting',globstring
        g=glob.glob(globstring)

        files=[]
        ras=[]
        decs=[]
        for d in g:
            hdu=fits.open(d)
            xsize=hdu[hdu_no].header['NAXIS1']
            ysize=hdu[hdu_no].header['NAXIS2']
            w=WCS(hdu[hdu_no].header)
            try:
                ra,dec=w.wcs_pix2world(xsize/2,ysize/2,0)
            except TypeError:
                ra,dec,_,_=w.wcs_pix2world(xsize/2,ysize/2,0,0,0)
            print ra,dec
            ras.append(ra)
            decs.append(dec)
            files.append(d)
            hdu.close()
        self.files=files
        self.ras=np.array(ras)
        self.decs=np.array(decs)
    def find(self,ra,dec):
        dist=(np.cos(self.decs*np.pi/180.0)*(self.ras-ra))**2.0 + (self.decs-dec)**2.0
        i=np.argmin(dist)
        return self.files[i]


if __name__=='__main__':
    t=Table.read(sys.argv[1])
    outfilename=sys.argv[1].replace('.fits','-list.txt')
    if not(os.path.isfile(outfilename)):
        startpoint=0
        outfile=open(outfilename,'w')
    else:
        startpoint=len(open(outfilename).readlines())
        outfile=open(outfilename,'a')
   
    imagedir=os.environ['IMAGEDIR']
    os.chdir(imagedir)

    lm=MapsList('mosaics/*-mosaic.fits')
    os.chdir('downloads')
    psm=MapsList('rings*.fits',hdu_no=1)
    wm=MapsList('*w1*fits')
    fm=MapsList('*%2B*.fits')

    for r in t[startpoint:]:
        lofarname=lm.find(r['RA'],r['DEC'])
        psnames=psm.find(r['RA'],r['DEC'])
        wisename=wm.find(r['RA'],r['DEC'])
        firstname=fm.find(r['RA'],r['DEC'])
        print >>outfile,r['Source_Name'],lofarname,psnames,wisename,firstname

    outfile.close()
