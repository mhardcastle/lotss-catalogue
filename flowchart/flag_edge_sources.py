#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
from astropy.table import Table,vstack,Column
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import sys
import os
import glob
#from subim import extract_subim
import multiprocessing as mp 
import itertools
import time
import matplotlib.pyplot as plt

IMAGEDIR='/data/lofar/mjh/hetdex_v4'

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import pyregion

def flatten(f,ra,dec,x,y,size,hduid=0,channel=0,freqaxis=3):
    """ 
    Flatten a fits file so that it becomes a 2D image. Return new header and
    data
    This version also makes a sub-image of specified size.
    """

    naxis=f[hduid].header['NAXIS']
    if naxis<2:
        raise OverlayException('Can\'t make map from this')

    #print f[hduid].data.shape
    ds=f[hduid].data.shape[-2:]
    by,bx=ds
    xmin=int(x-size)
    if xmin<0:
        xmin=0
    xmax=int(x+size)
    if xmax>bx:
        xmax=bx
    ymin=int(y-size)
    if ymin<0:
        ymin=0
    ymax=int(y+size)
    if ymax>by:
        ymax=by
    
    if ymax<=ymin or xmax<=xmin:
        # this can only happen if the required position is not on the map
        print 'Failed to make subimage!'
        print xmin,xmax,ymin,ymax
        return None

    w = WCS(f[hduid].header)
    wn=WCS(naxis=2)
    
    wn.wcs.crpix[0]=w.wcs.crpix[0]-xmin
    wn.wcs.crpix[1]=w.wcs.crpix[1]-ymin
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    try:
        wn.wcs.pc=w.wcs.pc[0:2,0:2]
    except AttributeError:
        pass # pc is not present
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]
    
    header = wn.to_header()
    header["NAXIS"]=2

    slice=[]
    for i in range(naxis,0,-1):
        if i==1:
            slice.append(np.s_[xmin:xmax+1])
        elif i==2:
            slice.append(np.s_[ymin:ymax+1])
        elif i==freqaxis:
            slice.append(channel)
        else:
            slice.append(0)
    #print slice

    hdu=fits.PrimaryHDU(f[hduid].data[slice],header)
    copy=('EQUINOX','EPOCH','BMAJ','BMIN','BPA')
    for k in copy:
        r=f[hduid].header.get(k)
        if r:
            hdu.header[k]=r
    if 'TAN' in hdu.header['CTYPE1']:
        hdu.header['LATPOLE']=f[hduid].header['CRVAL2']
    hdulist=fits.HDUList([hdu])
    return hdulist

def extract_subim(filename,ra,dec,size,hduid=0):
    #print 'Opening',filename
    orighdu=fits.open(filename)
    psize=int(size/orighdu[hduid].header['CDELT2'])
    ndims=orighdu[hduid].header['NAXIS']
    pvect=np.zeros((1,ndims))
    lwcs=WCS(orighdu[hduid].header)
    pvect[0][0]=ra
    pvect[0][1]=dec
    imc=lwcs.wcs_world2pix(pvect,0)
    x=imc[0][0]
    y=imc[0][1]
    hdu=flatten(orighdu,ra,dec,x,y,psize,hduid=hduid)
    return hdu


#def get_mosaic_name(name):
    #globst=IMAGEDIR+'/mosaics/'+name.rstrip()+'*'
    #g=glob.glob(globst)
    #if len(g)==1:
        #return g[0]
    #elif len(g)==0:
        #raise RuntimeError('No mosaic called '+name)
    #else:
        #raise RuntimeError('Mosaic name ambiguous')
    
def get_mosaic_name(name):
    if name == 'P21':
        name = 'P21-mosaic'
    g=glob.glob(IMAGEDIR+'/mosaics/'+name+'*')
    if len(g)>0:
        return g[0]
    else:
        raise RuntimeError('No mosaic called '+name)

    
    
def check_flagged(mosaicname,ra,dec,size):
    
    mosaicfile = get_mosaic_name(mosaicname)
    
    size = np.max((size, 0.001))
    lhdu = extract_subim(mosaicfile,ra,dec,size)
    
    
    if lhdu is None:
        print mosaicname,ra,dec,size
        return True
        
    
    data = lhdu[0].data
    
    flagged = np.any(np.isnan(data))
        
    
    return flagged


def check_flagged_region(mosaicname,ra,dec,maj,min,pa):
    
    mosaicfile = get_mosaic_name(mosaicname)
    
    maj = np.max((maj, 0.001))
    lhdu = extract_subim(mosaicfile,ra,dec,2*maj)
    
    
    if lhdu is None:
        print mosaicname,ra,dec,maj
        return True
        
    
    
    region = pyregion.parse('fk5;ellipse({ra},{dec},{a},{b},{pa})'.format(ra=ra,dec=dec,a=maj,b=min,pa=pa))
    
    data = lhdu[0].data
    region_mask = region.get_mask(lhdu[0])
    
    
    flagged = np.any(np.isnan(data[region_mask]))
    
        
    
    return flagged


def check_flagged_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return check_flagged(*a_b)
    
if __name__ == '__main__':
    lofarcat_file = '/data/lofar/wwilliams/hetdex/LoTSS-DR1-July21-2017/LOFAR_HBA_T1_DR1_catalog_v0.95_masked.srl.fits'
    lofarcat_flagged_file = '/data/lofar/wwilliams/hetdex/LoTSS-DR1-July21-2017/LOFAR_HBA_T1_DR1_catalog_v0.95_masked.srl.edgeflags.fits'
    
    lofarcat_all = Table.read(lofarcat_file)
    
    ncpu = 16
    pool = mp.Pool(processes=ncpu)
    
    #for some reason this stops midway?? do it in chunks
    
    Ns = ncpu
    N = ncpu*1200
    nanflags = np.zeros(len(lofarcat_all),dtype=bool)
    for i in range(Ns):
        i2 = np.min((len(lofarcat_all),(i+1)*N)) 
        print 'chunk {}: {}-{} '.format(i,i*N, i2)
        starttime = time.time()
        lofarcat = lofarcat_all[i*N:i2]
        results = pool.map(check_flagged_star, itertools.izip(lofarcat['Mosaic_ID'],lofarcat['RA'],lofarcat['DEC'],lofarcat['Maj']/3600.))
        nanflags[i*N:i2] = np.array(results)
        endtime = time.time()
    
        print 'Took {} seconds for {} sources' .format(endtime-starttime, len(lofarcat))
    
    
    
    # take a closer look at the ones that have been flagged as having a nanpixel within Maj
    # apply region mask and check for nan pixels inside that only
    nanflags2 = nanflags.copy()
    
    for ti in np.where(nanflags)[0]:
        tt = lofarcat_all[ti]
        newflag = check_flagged_region(tt['Mosaic_ID'],tt['RA'],tt['DEC'],tt['Maj']/3600.,tt['Min']/3600.,tt['PA'])
        print ti, newflag
        nanflags2[ti] = newflag
    
    
    lofarcat_all.add_column(Column(nanflags, 'Edge_flag1'))
    lofarcat_all.add_column(Column(nanflags2, 'Edge_flag2'))
    
    lofarcat_all.write(lofarcat_flagged_file,overwrite=True)
    
    
    # plot where the edges are
    plt.figure()
    plt.scatter(lofarcat_all['RA'],lofarcat_all['DEC'],c=1.*nanflags2,edgecolor='none')
    #plt.scatter(lofarcat['RA'][nanflags],lofarcat['DEC'][nanflags],c='y',edgecolor='none')
    plt.savefig('edgeflags.png')
    
    
    # check the ones highlighted in the issue
    checklist = ['ILTJ151524.99+543046.0',
                    'ILTJ150319.52+454944.5',
                    'ILTJ142450.82+502655.3',
                    'ILTJ142342.77+524831.0',
                    'ILTJ141236.02+512257.8',
                    'ILTJ124236.50+562844.6',
                    'ILTJ114646.72+562353.2']
    
    for n in checklist: 
        print n, n in lofarcat_all[nanflags2]['Source_Name']
        
    ## hmm, two not flagged - one is in a messy part of the sky, the other is near the edge but doesn't quite make the cut


    