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

    workflow=int(sys.argv[1])
    infile='/beegfs/lofar/deepfields/ELAIS_N1_LR/EN1_ML_RUN_fin_overlap_srl_workflow_fixed.fits'
    t=Table.read(infile)
    t=t[t['FLAG_WORKFLOW']==workflow]
    print len(t)
    t.write('workflow.fits',overwrite=True)
    ot=Table.read(infile)
    # large source table for neighbours
    lt=ot[(ot['DC_Maj']>8/3600.0)]

    # read in the big files that have all the data
    print 'Reading data...'
    allgals=Table.read('/beegfs/lofar/deepfields/ELAIS_N1_optical/catalogues/EN1_merged_pos_all.fits')
    gals=Table.read('/beegfs/lofar/deepfields/ELAIS_N1_optical/catalogues/EN1_merged_pos.fits')
    lofarfile=fits.open('/beegfs/lofar/deepfields/ELAIS-N1_LOFAR/image_full_ampphase_di_m.NS_shift.int.facetRestored.fits')
    spitzerfile=fits.open('/beegfs/lofar/deepfields/ELAIS_N1_optical/optical_images/sw2band/EL_EN1_sw2band.fits')
    ibandfile=fits.open('/beegfs/lofar/deepfields/ELAIS_N1_optical/optical_images/iband/EL_EN1_iband.fits')
    
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
        simage=sourcename+'_S.png'

        ra,dec=r['RA'],r['DEC']

        title=None
        maxsize=120
        minsize=60
        
        # resize the image to look for interesting neighbours
        iter=0
        while True:
            startra,startdec=ra,dec
            tcopy=lt
            tcopy['dist']=np.sqrt((np.cos(dec*np.pi/180.0)*(tcopy['RA']-ra))**2.0+(tcopy['DEC']-dec)**2.0)*3600.0
            tcopy=tcopy[tcopy['dist']<maxsize]
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
        
        if np.isnan(size):
            ra=r['RA']
            dec=r['DEC']
            size=minsize

        if size>maxsize:
            # revert just to original
            ra,dec=r['RA'],r['DEC']
            tcopy=Table(r)
            ra,dec,size=find_bbox(tcopy,scale=scale)

        if size>maxsize:
            size=maxsize
        if size<minsize:
            size=minsize
        size=(int(0.5+size/10))*10
        print 'size is',size

        size/=3600.0

        seps=separation(ra,dec,ot['RA'],ot['DEC'])
        ots=ot[seps<size]
        print ra,dec
        print ots['RA','DEC']

        ls=[]
        for nr in ots:
            if nr['Source_Name']==r['Source_Name']:
                ls.append('solid')
            else:
                ls.append('dashed')

        lhdu=extract_subim(lofarfile,ra,dec,size)
        shdu=extract_subim(spitzerfile,ra,dec,size)
        ihdu=extract_subim(ibandfile,ra,dec,size)

        pg=gals[(np.abs(gals['ra']-ra)<(size/np.cos(dec*np.pi/180.0))) & (np.abs(gals['dec']-dec)<size)]

        try:
            peak==r['Peak_flux']
        except:
            peak=None

        plist=[]
        if not np.isnan(r['lr_index_fin']):
            index=int(r['lr_index_fin'])
            marker_ra=allgals[index]['ra']
            marker_dec=allgals[index]['dec']
        else:
            marker_ra=None
            marker_dec=None
        if not os.path.isfile(iimage):
            show_overlay(lhdu,ihdu,ra,dec,size,firsthdu=None,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,save_name=iimage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,lw_ellipse=3,ellipse_style=ls,peak=peak,noisethresh=1,plotpos=pg,ppsize=350)
            plist.append(mogrify(iimage))
        if not os.path.isfile(simage):
            show_overlay(lhdu,shdu,ra,dec,size,firsthdu=None,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,save_name=simage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,lw_ellipse=3,ellipse_style=ls,peak=peak,plotpos=pg,ppsize=350)
            plist.append(mogrify(simage))

        #with open(manifestname,'w') as manifest:
        #    manifest.write('%i,%s,%s,%s,%s,%s,%f,%f,%f\n' % (i,iimage,ipimage,simage,spimage,sourcename,ra,dec,size*3600.0))

        for p in plist:
            p.wait()
            
        os.system('montage '+iimage+' '+simage+' -tile 2x1 -geometry 640x640+5+0 -background "#000000" '+sourcename+'_j.png')
