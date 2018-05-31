#!/usr/bin/python

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
import pickle
from separation import separation

class MapsList(object):
    def __init__(self,globstring,hdu_no=0,keephdu=True):
        warnings.simplefilter('ignore',category=AstropyWarning)
        print 'Ingesting',globstring
        g=glob.glob(globstring)

        files=[]
        ras=[]
        decs=[]
        sizes=[]
        wcs=[]
        hdus=[]
        for d in g:
            hdu=fits.open(d)
            xsize=hdu[hdu_no].header['NAXIS1']
            ysize=hdu[hdu_no].header['NAXIS2']
            w=WCS(hdu[hdu_no].header)
            try:
                ra,dec=w.wcs_pix2world(xsize/2,ysize/2,0)
            except TypeError:
                ra,dec,_,_=w.wcs_pix2world(xsize/2,ysize/2,0,0,0)
            ras.append(ra)
            decs.append(dec)
            files.append(d)
            wcs.append(w)
            sizes.append((xsize,ysize))
            if keephdu:
                hdus.append(hdu)
            else:
                hdu.close()
        self.files=files
        self.ras=np.array(ras)
        self.decs=np.array(decs)
        self.sizes=sizes
        self.wcs=wcs
        self.hdus=hdus
    def find(self,ra,dec,returnhdu=False):
        dist=separation(ra,dec,self.ras,self.decs)
        ranks=sorted(range(len(dist)),key=lambda i:dist[i])
        for r in ranks:
            w=self.wcs[r]
            xsize,ysize=self.sizes[r]
            try:
                x,y=w.wcs_world2pix(ra,dec,0)
            except TypeError:
                x,y,_,_=w.wcs_world2pix(ra,dec,0,0,0)
            if x>=0 and y>=0 and x<xsize and y<ysize:
                break
        else:
            raise RuntimeError('Cannot find suitable map')
        if returnhdu:
            return self.hdus[r]
        else:
            return self.files[r]

if __name__=='__main__':
    t=Table.read(sys.argv[1])
    outfilename=sys.argv[1].replace('.fits','-list.txt')
    outfile=open(outfilename,'w')
    
    imagedir=os.environ['IMAGEDIR']
    os.chdir(imagedir)
    if os.path.isfile('mapslist.pickle'):
        with open('mapslist.pickle') as f:
            lm,psm,wm,fm=pickle.load(f)
    else:
        lm=MapsList('mosaics/*-mosaic.fits')
        os.chdir('downloads')
        psm=MapsList('rings*.fits',hdu_no=1)
        wm=MapsList('*w1*fits')
        fm=MapsList('*%2B*.fits')
        os.chdir(imagedir)
        with open('mapslist.pickle','w') as f:
            pickle.dump([lm,psm,wm,fm],f)

    for r in t:
        lofarname=lm.find(r['RA'],r['DEC'])
        psnames=psm.find(r['RA'],r['DEC'])
        wisename=wm.find(r['RA'],r['DEC'])
        firstname=fm.find(r['RA'],r['DEC'])
        print >>outfile,r['Source_Name'],lofarname,psnames,wisename,firstname

    outfile.close()
