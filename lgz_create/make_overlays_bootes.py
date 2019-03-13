#!/usr/bin/python

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

scale=1.0 # scaling factor for region sizes

def mogrify(filename):
    command='mogrify -quality 90 -trim '+filename
    p = subprocess.Popen(command,shell=True)
    return p


if __name__=='__main__':

    tname=sys.argv[1]
    t=Table.read(tname)
    ot=Table.read('/beegfs/lofar/deepfields/Bootes_filter/Bootes_filtered_n.fits')
    # large source table for neighbours
    lt=ot[(ot['Total_flux']>0.003) & (ot['DC_Maj']>8/3600.0)]

    # read in the big files that have all the data
    print 'Reading data...'
    gals=Table.read('Bootes_MASTER_opt_spitzer_merged.fits')
    gals=gals[gals['FLAG_DEEP']==1]
    gals['ALPHA_J2000'].name='ra'
    gals['DELTA_J2000'].name='dec'
    lofarfile=fits.open('/beegfs/lofar/deepfields/Bootes_LOFAR/image_full_ampphase_di_m.NS_shift.int.facetRestored.blanked.scaled.fits')
    spitzerfile=fits.open('/beegfs/lofar/deepfields/Bootes_optical/SDWFS/I2_bootes.v32.fits')
    iband=ReadMaps('/beegfs/lofar/deepfields/Bootes_optical/NDWFS/I/*_con.fits')
    
    start=int(sys.argv[2])
    try:
        end=int(sys.argv[3])+1
    except:
        end=start+1
    if end>len(t):
        end=len(t)
        
    for i in range(start,end):

        r=t[i]
        print r
        sourcename=r['Source_Name']

        iimage=sourcename+'_I.png'
        ipimage=sourcename+'_Ip.png'
        simage=sourcename+'_S.png'
        manifestname=sourcename+'-manifest.txt'

        '''
        if os.path.isfile(manifestname):
            print 'Selected output file exists already'
            continue
        '''
        ra,dec=r['RA'],r['DEC']

        marker_ra=None
        marker_dec=None
        title=None

        # resize the image to look for interesting neighbours
        iter=0
        while True:
            startra,startdec=ra,dec
            tcopy=lt
            tcopy['dist']=np.sqrt((np.cos(dec*np.pi/180.0)*(tcopy['RA']-ra))**2.0+(tcopy['DEC']-dec)**2.0)*3600.0
            tcopy=tcopy[tcopy['dist']<90]
            print 'Iter',iter,'found',len(tcopy),'neighbours'

            # make sure the original source is in there
            for nr in tcopy:
                if sourcename==nr['Source_Name']:
                    break
            else:
                if 'Maj' in r.columns:
                    tcopy=vstack((tcopy,r))

            ra=np.mean(tcopy['RA'])
            dec=np.mean(tcopy['DEC'])

            if startra==ra and startdec==dec:
                break
            iter+=1
            print iter,len(tcopy)
            if iter==10:
                break

            # now find the bounding box of the resulting collection
        ra,dec,size=find_bbox(tcopy,scale=scale)
        size/=1.5 # scale down
        
        if np.isnan(size):
            ra=r['RA']
            dec=r['DEC']
            size=30

        if size>90.0:
            # revert just to original
            ra,dec=r['RA'],r['DEC']
            tcopy=Table(r)
            ra,dec,size=find_bbox(tcopy,scale=scale)

        if size>90:
            size=90.0
        if size<30:
            size=30.0
        size=(int(0.5+size/10))*10
        print 'size is',size

        size/=3600.0

        seps=separation(ra,dec,ot['RA'],ot['DEC'])
        ots=ot[seps<size*2]
        print ra,dec
        print ots['RA','DEC']

        cols=[]
        for nr in ots:
            if nr['Source_Name']==r['Source_Name']:
                # this is our target
                cols.append('red')
            else:
                cols.append('green')
        # here we use Montage to make a regridded Spitzer image so that
        # the three images look the same
        
        lhdu=extract_subim(lofarfile,ra,dec,size)
        shdu=extract_subim(spitzerfile,ra,dec,size*1.5)
        imain=iband.find_pos(ra,dec,hdu=True)
        ihdu=extract_subim(imain,ra,dec,size)
        ihdu.writeto(sourcename+'_i.fits',overwrite=True)
        shdu.writeto(sourcename+'_s.fits',overwrite=True)
        montage_wrapper.mGetHdr(sourcename+'_i.fits',sourcename+'_i.hdr') 
        montage_wrapper.mProject(sourcename+'_s.fits',sourcename+'_so.fits',sourcename+'_i.hdr')
        ihdu[0].data=np.where(ihdu[0].data>49999,np.nan,ihdu[0].data)
        shdu=fits.open(sourcename+'_so.fits')

        pg=gals[(np.abs(gals['ra']-ra)<(1.2*size/np.cos(dec*np.pi/180.0))) & (np.abs(gals['dec']-dec)<size)]

        try:
            peak==r['Peak_flux']
        except:
            peak=None

        print ots,cols

        plist=[]
        if not os.path.isfile(iimage):
            show_overlay(lhdu,ihdu,ra,dec,size,firsthdu=None,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,save_name=iimage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,ellipse_color=cols,peak=peak,noisethresh=1)
            plist.append(mogrify(iimage))
        if not os.path.isfile(ipimage):
            show_overlay(lhdu,ihdu,ra,dec,size,firsthdu=None,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,save_name=ipimage,show_lofar=False,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,ellipse_color=cols,peak=peak,noisethresh=1,plotpos=pg,ppsize=350)
            plist.append(mogrify(ipimage))
        if not os.path.isfile(simage):
            show_overlay(lhdu,shdu,ra,dec,size,firsthdu=None,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,save_name=simage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,ellipse_color=cols,peak=peak)
            plist.append(mogrify(simage))

        with open(manifestname,'w') as manifest:
            manifest.write('%i,%s,%s,%s,%s,%f,%f,%f\n' % (i,iimage,ipimage,simage,sourcename,ra,dec,size*3600.0))

        for p in plist:
            p.wait()
            
