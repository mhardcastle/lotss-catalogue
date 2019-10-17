import astropy
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

def flatten(f,ra,dec,x,y,size,hduid=0,channel=0,freqaxis=3,verbose=True):
    """ 
    Flatten a fits file so that it becomes a 2D image. Return new header and
    data
    This version also makes a sub-image of specified size.
    """

    naxis=f[hduid].header['NAXIS']
    if naxis<2:
        raise RuntimeError('Can\'t make map from this')

    if verbose:
        print 'Input image shape is',f[hduid].data.shape
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
        print xmin,xmax,ymin,ymax
        raise RuntimeError('Failed to make subimage!')

    w = WCS(f[hduid].header)
    wn=WCS(naxis=2)
    
    wn.wcs.crpix[0]=w.wcs.crpix[0]-xmin
    wn.wcs.crpix[1]=w.wcs.crpix[1]-ymin
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    '''
    try:
        wn.wcs.pc=w.wcs.pc[0:2,0:2]
    except AttributeError:
        pass # pc is not present
    try:
        wn.wcs.cd=w.wcs.cd[0:2,0:2]
    except AttributeError:
        pass # cd is not present
    '''
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]
    
    header = wn.to_header()
    header["NAXIS"]=2

    slice=[]
    for i in range(naxis,0,-1):
        if i==1:
            slice.append(np.s_[xmin:xmax])
        elif i==2:
            slice.append(np.s_[ymin:ymax])
        elif i==freqaxis:
            slice.append(channel)
        else:
            slice.append(0)
    if verbose:
        print slice

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

def extract_subim(filename,ra,dec,size,hduid=0,verbose=True):
    """Extract a sub-image and return an HDU.
    filename: the input FITS file
    ra, dec: the position in degrees
    size: the half-size in degrees
    hduid: the element of the original HDU to use (default 0)
    verbose: print diagnostics (default True)
    """
    if isinstance(filename,astropy.io.fits.hdu.hdulist.HDUList):
        orighdu=filename
    else:
        if verbose:
            print 'Opening',filename
        orighdu=fits.open(filename)
    if 'CDELT2' in orighdu[hduid].header:
        delt=orighdu[hduid].header['CDELT2']
    else:
        # assuming no rotation here
        delt=orighdu[hduid].header['CD2_2']
        
    psize=int(size/delt)
    
    ndims=orighdu[hduid].header['NAXIS']
    pvect=np.zeros((1,ndims))
    lwcs=WCS(orighdu[hduid].header)
    pvect[0][0]=ra
    pvect[0][1]=dec
    imc=lwcs.wcs_world2pix(pvect,0)
    x=imc[0][0]
    y=imc[0][1]
    if verbose:
        print 'Extracting sub-image'
    hdu=flatten(orighdu,ra,dec,x,y,psize,hduid=hduid,verbose=verbose)
    '''
    del(hdu[hduid].header['PC1_1'])
    del(hdu[hduid].header['PC2_2'])
    '''
    hdu[hduid].header['CDELT2']=delt
    hdu[hduid].header['CDELT1']=-delt
    return hdu
