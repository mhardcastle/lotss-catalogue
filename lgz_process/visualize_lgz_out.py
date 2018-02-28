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

scale=3600.0 # scaling factor for region sizes

if __name__=='__main__':

    imagedir=os.environ['IMAGEDIR']
    tname=sys.argv[1]
    lname=sys.argv[1].replace('.fits','-list.txt')
    t=Table.read(tname)
    # Annotated PyBDSF table
    ot=Table.read('/data/lofar/wwilliams/hetdex/LoTSS-DR1-July21-2017/LOFAR_HBA_T1_DR1_merge_ID_v0.8.comp.fits')
    # ID catalogue so we can check for a source with existing ML or other ID
    it=Table.read('/data/lofar/mjh/hetdex_v4/lgz_v2_output/HETDEX-LGZ-comps-lgzv2.fits')
    for n in it.colnames:
        if 'Comp_' in n:
            it[n].name=n.replace('Comp_','')
    # read lists
    lines=[l.rstrip().split() for l in open(lname).readlines()]
    lofarmaps=[l[1] for l in lines]
    psmaps=[l[2] for l in lines]
    wisemaps=[l[3] for l in lines]
    firstmaps=[l[4] for l in lines]

#    lofarmaps=[l.rstrip().split()[1] for l in open('lofar-maps.txt').readlines()]
#    psmaps=[l.rstrip().split()[1] for l in open('panstarrs-list.txt').readlines()]
#    firstmaps=[l.rstrip().split()[1] for l in open('first-list.txt').readlines()]
#    wisemaps=[l.rstrip().split()[2] for l in open('panstarrs-list.txt').readlines()]

    i=int(sys.argv[2])
    r=t[i]
    print r
    sourcename=r['Source_Name']
    sourcename=sourcename.rstrip()

    psimage=sourcename+'_PS.png'
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


    ra,dec=r['RA'],r['DEC']

    marker_ra=r['optRA']
    marker_dec=r['optDec']

    title=sourcename

    # look in component catalogue for associated comps
    mask=(it['Source_Name']==sourcename)
    tcopy=it[mask]
    assert(len(tcopy)>0)
    ra,dec,size=find_bbox(tcopy)
    size*=1.5

    
    if np.isnan(size):
        ra=r['RA']
        dec=r['DEC']
        size=60

    if size>300.0:
        # revert just to original
        ra,dec=r['RA'],r['DEC']
        tcopy=Table(r)
        ra,dec,size=find_bbox(tcopy)

    if size>300:
        size=300.0
    if size<60:
        size=60.0
    size=(int(0.5+size/10))*10
    print 'size is',size

    size/=3600.0

    seps=separation(ra,dec,ot['RA'],ot['DEC'])
    ots=ot[seps<size*2]
    print ra,dec
    print ots['RA','DEC']
    ots=ots[ots['New_Source_Name']!=""] # removes artefacts
    cols=[]
    for nr in ots:
        name=nr['Source_Name']
        for c in tcopy:
            if c['Name']==name:
                cols.append('green')
                break
        else:
            cols.append('red')
    
    pshdu=extract_subim(imagedir+'/downloads/'+psmaps[i],ra,dec,size,hduid=1)
    lhdu=extract_subim(lofarfile,ra,dec,size)
    firsthdu=extract_subim(imagedir+'/downloads/'+firstmaps[i],ra,dec,size)
    try:
        peak==r['Peak_flux']/1000.0
    except:
        peak=None

    print ots,cols
    
    show_overlay(lhdu,pshdu,ra,dec,size,firsthdu=firsthdu,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,save_name=psimage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,ellipse_color=cols,peak=peak)
