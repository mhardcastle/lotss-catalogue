#!/usr/bin/python

# version that'll use the source and component catalogues already generated

from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from astropy.table import Table,vstack
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import sys
import os
import glob
from subim import extract_subim
from separation import separation
from image_utils import find_bbox,get_mosaic_name
from overlay import show_overlay
from read_maps import ReadMaps
import subprocess
import montage_wrapper

scale=3600.0 # scaling factor for region sizes

def mogrify(filename):
    command='mogrify -quality 90 -trim '+filename
    p = subprocess.Popen(command,shell=True)
    return p


if __name__=='__main__':

    tname=sys.argv[1]
    t=Table.read(tname)
    lname=tname.replace('.fits','-list.txt')
    lines=[l.rstrip().split() for l in open(lname).readlines()]
    lofarmaps=[l[1] for l in lines]
    wisemaps=[l[2] for l in lines]
 
    version='v0.7'
    wd='/beegfs/lofar/mjh/rgz/Spring'
    
    dir=os.getcwd()
    print('Reading data...')
    

    st=Table.read(wd+'/sources-'+version+'.fits')
    ot=Table.read(wd+'/components-'+version+'.fits')

    # read in the big files that have all the data
    #gals=Table.read(wd+'/optical.fits')
    
    #gals['RA'].name='ra'
    #gals['DEC'].name='dec'
    # large source table for neighbours
    lt=ot[(ot['DC_Maj']>8/3600.0)]

    start=int(sys.argv[2])
    try:
        end=int(sys.argv[3])+1
    except:
        end=start+1
    if end>len(t):
        end=len(t)
        
    for i in range(start,end):

        r=t[i]
        print(r)
        sourcename=r['Source_Name']
        # find components of this source
        cfilter=ot['Parent_Source']==sourcename
        ctable=ot[cfilter]
        if len(ctable)==0:
            print('Source has no components!')
            continue
        
        simage=sourcename+'_S.png'
        spimage=sourcename+'_Sp.png'
        wisefile=os.environ['IMAGEDIR']+'/downloads/'+wisemaps[i]
        lofarfile=os.environ['IMAGEDIR']+'/'+lofarmaps[i]
        
        if os.path.isfile(simage):
            print('Selected output file exists already')
            continue

        ra,dec=r['RA'],r['DEC']

        marker_ra=None
        marker_dec=None
        title=None
        maxsize=300
        minsize=60
        
        # resize the image to look for interesting neighbours
        iter=0
        while True:
            startra,startdec=ra,dec
            tcopy=lt
            tcopy['dist']=np.sqrt((np.cos(dec*np.pi/180.0)*(tcopy['RA']-ra))**2.0+(tcopy['DEC']-dec)**2.0)*3600.0
            tcopy=tcopy[tcopy['dist']<maxsize]
            print('Iter',iter,'found',len(tcopy),'neighbours')

            # make sure the original source is in there
            tcopy=vstack((tcopy,ctable))
            print(tcopy)
            ra=np.mean(tcopy['RA'])
            dec=np.mean(tcopy['DEC'])

            if startra==ra and startdec==dec:
                print('converged in ra and dec')
                
                break
            iter+=1
            print(iter,len(tcopy))
            if iter==10:
                break

        # now find the bounding box of the resulting collection
        ra,dec,size=find_bbox(tcopy,scale=scale)
        
        if np.isnan(size):
            ra=r['RA']
            dec=r['DEC']
            size=minsize

        if size>maxsize:
            # revert just to original
            ra,dec=r['RA'],r['DEC']
            tcopy=ctable
            ra,dec,size=find_bbox(tcopy,scale=scale)

        if size>maxsize:
            size=maxsize
        if size<minsize:
            size=minsize
        size=(int(0.5+size/10))*10
        print('size is',size)

        size/=3600.0

        seps=separation(ra,dec,ot['RA'],ot['DEC'])
        ots=ot[seps<size]
        print(ra,dec)
        print(ots['RA','DEC'])

        ls=[]
        cs=[]
        for nr in ots:
            if nr['Parent_Source']==sourcename:
                ls.append('solid')
                cs.append('green')
            else:
                ls.append('dashed')
                cs.append('red')

        # here we use Montage to make a regridded Spitzer image so that
        # the three images look the same
        
        lhdu=extract_subim(lofarfile,ra,dec,size)
        try:
            shdu=extract_subim(wisefile,ra,dec,size)
        except RuntimeError as e:
            print('*** extract_subim failed (%s) ***' % str(e))
            continue
        
        #ihdu.writeto(sourcename+'_i.fits',overwrite=True)
        #shdu.writeto(sourcename+'_s.fits',overwrite=True)
        #montage_wrapper.mGetHdr(sourcename+'_i.fits',sourcename+'_i.hdr') 
        #montage_wrapper.mProject(sourcename+'_s.fits',sourcename+'_so.fits',sourcename+'_i.hdr')
        #ihdu[0].data=np.where(ihdu[0].data>49999,np.nan,ihdu[0].data)
        #shdu=fits.open(sourcename+'_so.fits')

        #pg=gals[(np.abs(gals['ra']-ra)<(size/np.cos(dec*np.pi/180.0))) & (np.abs(gals['dec']-dec)<size)]

        try:
            peak==r['Peak_flux']
        except:
            peak=None

        if not np.isnan(r['optRA']):
            marker_ra=r['optRA']
            marker_dec=r['optDec']
        else:
            marker_ra=None
            marker_dec=None

        print('Marker is going at',marker_ra,marker_dec)
        plist=[]
        error=False
        if not os.path.isfile(simage):
            try:
                show_overlay(lhdu,shdu,ra,dec,size,firsthdu=None,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,save_name=simage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='magenta',title=title,lw_ellipse=3,ellipse_style=ls,ellipse_color=cs,peak=peak,ppsize=350)
            except (Exception, RuntimeError) as e:
                error=True
            if error:
                print('*** image build failed! (%s) ***' % str(e))
                continue
            else:
                os.system('mogrify -quality 90 -trim '+sourcename+'*.png')
