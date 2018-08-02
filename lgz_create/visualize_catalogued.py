#!/usr/bin/python

# visualize sources in the final LGZ catalogue

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

scale=3600.0 # scaling factor for region sizes

if __name__=='__main__':

    imagedir=os.environ['IMAGEDIR']
    tname=sys.argv[1]
    lname=sys.argv[1].replace('.fits','-list.txt')
    t=Table.read(tname)
    # Annotated PyBDSF table
    ot=Table.read('/beegfs/general/lofar/LOFAR_HBA_T1_DR1_merge_ID_v1.0.comp.fits')
    # large source table for neighbours
    lt=ot[(ot['Total_flux']>3) & (ot['Maj']>8)]
    # ID catalogue so we can check for a source with existing ML or other ID
    it=Table.read('/beegfs/general/lofar/LOFAR_HBA_T1_DR1_merge_ID_optical_v1.0.fits')
    # read lists
    lines=[l.rstrip().split() for l in open(lname).readlines()]
    lofarmaps=[l[1] for l in lines]
    psmaps=[l[2] for l in lines]
    wisemaps=[l[3] for l in lines]
    firstmaps=[l[4] for l in lines]
    try:
        vlassmaps=[l[5] for l in lines]
    except:
        vlassmaps=['None']*len(firstmaps)

#    lofarmaps=[l.rstrip().split()[1] for l in open('lofar-maps.txt').readlines()]
#    psmaps=[l.rstrip().split()[1] for l in open('panstarrs-list.txt').readlines()]
#    firstmaps=[l.rstrip().split()[1] for l in open('first-list.txt').readlines()]
#    wisemaps=[l.rstrip().split()[2] for l in open('panstarrs-list.txt').readlines()]

    begin=int(sys.argv[2])
    try:
        end=int(sys.argv[3])
    except:
        end=begin
    if end>len(t):
        end=len(t)
        
    for i in range(begin,end):
    
        r=t[i]
        print r
        try:
            sourcename=r['Source_Name']
        except KeyError:
            # try to fix up old-format files
            t['Source_id'].name='Source_Name'
            sourcename=r['Source_Name']
        sourcename=sourcename.rstrip()

        psimage=sourcename+'_PS.png'
        pspimage=sourcename+'_PSp.png'
        wiseimage=sourcename+'_W.png'
        print lofarmaps[i],psmaps[i],wisemaps[i]
        try:
            mosaic=r['Mosaic_ID']
        except:
            mosaic=None
        if mosaic is not None:
            try:
                lofarfile=get_mosaic_name(mosaic)
            except RuntimeError:
                print 'Could not get mosaic ID from',r['Mosaic_ID']
                mosaic=None
        if mosaic is None:
            lofarfile=os.environ['IMAGEDIR']+'/'+lofarmaps[i]
        if os.path.isdir(lofarfile):
            lofarfile+='/mosaic.fits'
        print lofarfile

        from overlay import show_overlay

        ra,dec=r['RA'],r['DEC']

        try:
            marker_ra=r['ID_ra']
            marker_dec=r['ID_dec']
        except:
            marker_ra=None
            marker_dec=None

        title=None

        tmask=(ot['Source_Name']==sourcename)
        if np.sum(tmask)>0:
            tcopy=ot[tmask]
            ra,dec,size=find_bbox(tcopy)
            size*=1.2
        else:
            size=np.nan

        if np.isnan(size):
            ra=r['RA']
            dec=r['DEC']
            size=60

        if size<60:
            size=60.0
        size=(int(0.5+size/10))*10
        print 'size is',size

        size/=3600.0

        seps=separation(ra,dec,ot['RA'],ot['DEC'])
        ots=ot[seps<size*2]
        print ra,dec
        print ots['RA','DEC']
        ots=ots[ots['Source_Name']!=""] # removes artefacts
        cols=[]
        peak=None
        for nr in ots:
            if nr['Source_Name']==r['Source_Name']:
                # this is our target
                cols.append('green')
                if peak is None or nr['Peak_flux']<peak:
                    peak=nr['Peak_flux']
            else:
                cols.append('red')

        if 'NO_OVERLAY' in os.environ:
            ots=None

                
        pshdu=extract_subim(imagedir+'/downloads/'+psmaps[i],ra,dec,size,hduid=1)
        lhdu=extract_subim(lofarfile,ra,dec,size)
        firsthdu=extract_subim(imagedir+'/downloads/'+firstmaps[i],ra,dec,size)
        if vlassmaps[i]!='None':
            vlasshdu=extract_subim(imagedir+'/downloads/'+vlassmaps[i],ra,dec,size)
        else:
            vlasshdu=None

        if peak is not None:
            peak/=1000
        else:
            try:
                peak==r['Peak_flux']/1000.0
            except:
                peak=None

        print ots,cols
        try:
            plotpos=os.environ['PLOTPOS']
        except:
            plotpos=None
        if plotpos is not None:
            print 'PLOTPOS is',plotpos
            bits=plotpos.split(',')
            tp=Table.read(bits[0])
            try:
                tp['RA'].name='ra'
                tp['DEC'].name='dec'
            except:
                pass
            tpp=tp[(np.abs(tp['ra']-ra)<size) & (np.abs(tp['dec']-dec)<size)]
            print '*** Filtered to',len(tpp)
            shape='x'
            color='white'
            if len(bits)>1:
                shape=bits[1]
            if len(bits)>2:
                color=bits[2]
            plotpos=[(tpp,shape)]
            print '*** PLOTPOS is ',plotpos

        show_overlay(lhdu,pshdu,ra,dec,size,firsthdu=firsthdu,vlasshdu=vlasshdu,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,save_name=psimage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='magenta',title=title,ellipse_color=cols,peak=peak,plotpos=plotpos)
        whdu=extract_subim(imagedir+'/downloads/'+wisemaps[i],ra,dec,size)
        show_overlay(lhdu,whdu,ra,dec,size,firsthdu=firsthdu,vlasshdu=vlasshdu,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,save_name=wiseimage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='magenta',title=title,noisethresh=0,ellipse_color=cols,peak=peak,plotpos=plotpos)
