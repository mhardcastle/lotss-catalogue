#!/usr/bin/python

# download required image files (WISE, PanSTARRS, FIRST) for a set of
# input objects in a FITS catalogue. Outputs an image file and populates the
# $IMAGEDIR/downloads directory

import requests
from astropy.io import fits
from astropy.table import Table
import astropy.coordinates as coord
import astropy.units as u
from astroquery.ibe import Ibe
import shutil
import os.path
import sys
import re
import numpy as np
from lxml import html
import glob
import os
from time import sleep

def download_file(url,outname):
    if os.path.isfile(outname):
        print 'File',outname,'already exists, skipping'
    else:
        print 'Downloading',outname
        while True:
            try:
                response = requests.get(url, stream=True,verify=False,timeout=120)
                if response.status_code!=200:
                    print 'Warning, HTML status code',response.status_code
                    if response.status_code>=500 and response.status_code<600:
                        raise RuntimeError('Retry!')
            except requests.exceptions.ConnectionError:
                print 'Connection error! sleeping 60 seconds before retry...'
                sleep(60)
            except RuntimeError:
                print 'Transient error reported, retrying after sleep'
                sleep(60)
            except requests.exceptions.Timeout:
                print 'Timeout: retrying download'
            else:
                break
        with open(outname, 'wb') as out_file:
            shutil.copyfileobj(response.raw, out_file)
        del response

def get_panstarrs(ra,dec,psband):
    page=requests.get('http://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra=%f&dec=%f' % (ra,dec),verify=False)
    print page.status_code
    print page.headers['content-type']
    lines=page.text.split('\n')
    downloads=[]
    for l in lines[1:]:
        bits=l.split()
        if bits[4] in psband or psband=='*':
            download_file('http://ps1images.stsci.edu/'+bits[7],bits[8])
            downloads.append(bits[8])
    return downloads

def get_wise(ra,dec,band):
    mission='wise'
    dataset='allwise'
    table='p3am_cdd'
    successful=False
    while not successful:
        try:
            results=Ibe.query_region(coord.SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs'),mission=mission,dataset=dataset,table=table)
            successful=True
        except requests.exceptions.ConnectionError:
            print 'Connection failed, retrying'
            sleep(10)
    
    url = 'http://irsa.ipac.caltech.edu/ibe/data/'+mission+'/'+dataset+'/'+table+'/'
    params = { 'coadd_id': results[results['band']==band]['coadd_id'][0],
           'band': band }
    params['coaddgrp'] = params['coadd_id'][:2]
    params['coadd_ra'] = params['coadd_id'][:4]
    path = str.format('{coaddgrp:s}/{coadd_ra:s}/{coadd_id:s}/{coadd_id:s}-w{band:1d}-int-3.fits',**params)
    outname=path.split('/')[-1]
    download_file(url+path,outname)
    return outname

def get_first(ra,dec):
    url="http://archive.stsci.edu/"
    page=requests.get(url+"vlafirst/search.php?RA=%.7f&DEC=%.6f&Radius=30.0&action=Search" % (ra,dec),verify=False)
    print page.status_code

    tree=html.fromstring(page.text)
    table=tree.xpath('//tbody')
    links=[]
    dists=[]
    for row in table[0].getchildren():
        td=row.getchildren()
        links.append(td[0].getchildren()[0].attrib['href'])
        dists.append(float(td[8].text))

    index=np.argmin(dists)
    path=links[index]
   
    outname=path.split('/')[-1]
    download_file(url+path,outname)
    return outname

def get_legacy(ra,dec,size=1000,pixscale=0.454,bands='r',ftype='fits'):
    url = "http://legacysurvey.org/viewer/{}-cutout?ra={}&dec={}&size={}&layer=mzls+bass-dr6&pixscale={}&bands={}".format(ftype,ra,dec,size,pixscale,bands)
    outname='legacy-%s-%f-%f.fits' % (bands,ra,dec)
    download_file(url,outname)
    return outname
    
class LofarMaps(object):
    def __init__(self,stay_in_imagedir=False):
        imagedir=os.environ['IMAGEDIR']
        wd=os.getcwd()
        os.chdir(imagedir)
        # read the LOFAR map positions
        g=glob.glob('mosaics/P*')

        files=[]
        ras=[]
        decs=[]
        for d in g:
            if '.fits' in d:
                file=d
            else:
                file=d+'/mosaic.fits'
                if not os.path.isfile(file):
                    continue
            hdu=fits.open(file)
            ras.append(hdu[0].header['CRVAL1'])
            decs.append(hdu[0].header['CRVAL2'])
            files.append(file)
        self.files=files
        self.ras=np.array(ras)
        self.decs=np.array(decs)
        if not stay_in_imagedir:
            os.chdir(wd)
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

    lm=LofarMaps(stay_in_imagedir=True)
    
    # now work in downloads dir
    os.chdir('downloads')

    for r in t[startpoint:]:
        lofarname=lm.find(r['RA'],r['DEC'])
        psnames=get_panstarrs(r['RA'],r['DEC'],'i')
        wisename=get_wise(r['RA'],r['DEC'],1)
        firstname=get_first(r['RA'],r['DEC'])
        print >>outfile,r['Source_Name'],lofarname,psnames[0],wisename,firstname

    outfile.close()
