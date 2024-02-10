# VLASS downloader

from __future__ import print_function
import requests
import re
from separation import separation
from download_file import download_file
from reproject import reproject_interp
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from copy import deepcopy
import os

def flatten(f):
    """ Flatten a fits file so that it becomes a 2D image. Return new header and data. From auxcodes.py """

    naxis=f[0].header['NAXIS']
    if naxis<2:
        raise RuntimeError('Can\'t make map from this')
    if naxis==2:
        f[0].header["WCSAXES"]=2
        return fits.PrimaryHDU(header=f[0].header,data=f[0].data)

    w = WCS(f[0].header)
    wn=WCS(naxis=2)
    
    wn.wcs.crpix[0]=w.wcs.crpix[0]
    wn.wcs.crpix[1]=w.wcs.crpix[1]
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]
    
    header = wn.to_header()
    header["NAXIS"]=2
    header["WCSAXES"]=2
    copy=('EQUINOX','EPOCH','BMAJ', 'BMIN', 'BPA', 'RESTFRQ', 'TELESCOP', 'OBSERVER')
    for k in copy:
        r=f[0].header.get(k)
        if r is not None:
            header[k]=r

    dslice=[]
    for i in range(naxis,0,-1):
        if i<=2:
            dslice.append(np.s_[:],)
        else:
            dslice.append(0)
        
    hdu = fits.PrimaryHDU(header=header,data=f[0].data[tuple(dslice)])
    return hdu

def get_vlass_tiles(ra,dec):
    # Find a position on the sky and see which relevant tiles have been imaged
    page=requests.get('https://archive-new.nrao.edu/vlass/VLASS_dyn_summary.php')
    print(page.status_code)
    print(page.headers['content-type'])
    lines=page.text.split('\n')
    results=[]
    
    for l in lines[3:]:
        if not l: continue
        bits=l.split()
        dec_min=float(bits[1])
        dec_max=float(bits[2])
        ra_min=float(bits[3])*15
        ra_max=float(bits[4])*15
        if ra>=ra_min and ra<=ra_max and dec>dec_min and dec<dec_max and 'imaged' in l:
            results.append((bits[0],bits[5]))
    return results

def get_vlass(ra,dec):
    files=[]
    tiles=get_vlass_tiles(ra,dec)
    if len(tiles)==0:
        raise RuntimeError('No VLASS coverage of this position')
    for tile,epoch in tiles:
        print(tile,epoch)
        if epoch.startswith('VLASS1'):
            epoch+='v2'
        urlbase='https://archive-new.nrao.edu/vlass/quicklook/'+epoch+'/'+tile
        print('Checking',urlbase)
        page=requests.get(urlbase)
        lines=page.text.split('\n')
        wanted_fields=[]
        for l in lines:
            if '<td>' in l: # a table line
                dirname=re.sub('/.*','/',re.sub('.*VLASS','VLASS',l))
                if 'J' not in dirname: continue # upward link
                pos=re.sub('\..*','',re.sub('.*J','J',dirname))
                ra_h=int(pos[1:3])
                ra_m=int(pos[3:5])
                ra_s=int(pos[5:7])
                decsign=pos[7]
                dec_d=int(pos[8:10])
                dec_m=int(pos[10:12])
                dec_s=int(pos[12:14])
                i_ra=15.0*(ra_h+ra_m/60.0+ra_s/3600.0)
                i_dec=dec_d+dec_m/60.0+dec_s/3600.0
                if decsign=='-': i_dec=-i_dec
                sep=separation(ra,dec,i_ra,i_dec)
                #print(pos,ra_h,ra_m,ra_s,decsign,dec_d,dec_m,dec_s,i_ra,i_dec,separation(ra,dec,i_ra,i_dec))
                if sep<1.0:
                    wanted_fields.append(dirname)

        for field in wanted_fields:
            filename=field.replace('/','.I.iter1.image.pbcor.tt0.subim.fits')
            download_file(urlbase+'/'+field+filename,filename)
            files.append(filename)

    return files

def mosaic_vlass(ra,dec,size=1000,outname=None,overwrite=False):
    if outname is None:
        outname='VLASS-mosaic-%.2f_%.2f.fits' % (ra,dec)
    if os.path.isfile(outname) and not overwrite:
        print('Output file %s  exists, skipping!' % outname)
        return False
    files=get_vlass(ra,dec)
    if len(files)==0:
        raise RuntimeError('No VLASS files for this position!')
    hdus=[]
    for f in files:
        hdu=fits.open(f)
        hdus.append(flatten(hdu))

    header=deepcopy(hdus[0].header)
    header['NAXIS1']=size
    header['NAXIS2']=size
    header['CRPIX1']=size/2
    header['CRPIX2']=size/2
    header['CRVAL1']=ra
    header['CRVAL2']=dec

    isum=np.zeros([size,size])
    wsum=np.zeros_like(isum) # will be 0,1,2 etc...
    mask=np.zeros_like(isum,dtype=np.bool) # True if there is at least one non-NaN contribution

    for i in range(len(hdus)):
        print('image',i,'(',files[i],')')
        r=reproject_interp(hdus[i],header,return_footprint=False)
        print(np.sum(~np.isnan(r)),'non-nan pix')
        w=np.where(np.isnan(r),0,1)
        mask|=~np.isnan(r)
        r[np.isnan(r)]=0
        isum+=r
        wsum+=w

    isum/=wsum
    isum[~mask]=np.nan # blank
    hdu = fits.PrimaryHDU(header=header,data=isum)
    if outname is None:
        outname='VLASS-mosaic-%.2f_%.2f.fits' % (ra,dec)
    hdu.writeto(outname,overwrite=True)
    return True
    
if __name__=='__main__':
    print('Files downloaded are',get_vlass(208.07452,31.44625))
    
