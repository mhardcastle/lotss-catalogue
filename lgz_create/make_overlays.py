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
from image_utils import find_bbox,get_mosaic_name

scale=3600.0 # scaling factor for region sizes


if __name__=='__main__':

    imagedir=os.environ['IMAGEDIR']
    tname=sys.argv[1]
    lname=sys.argv[1].replace('.fits','-list.txt')
    t=Table.read(tname)
    # all LOFAR sources
    ot=Table.read(imagedir+'/LOFAR_HBA_T1_DR1_catalog_v0.2.srl.fits')
    # large source table
    lt=ot[(ot['Total_flux']>10) & (ot['Maj']>8)]

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
    manifestname=sourcename+'-manifest.txt'
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

    if os.path.isfile(manifestname):
        print 'Selected output file exists already'
        sys.exit(0)

    from overlay import show_overlay

    ra,dec=r['RA'],r['DEC']

    try:
        marker_ra=r['ra']
        marker_dec=r['dec']
    except:
        marker_ra=None
        marker_dec=None

    try:
        title='LR = %f' % r['lr_pc_7th']
    except:
        title=None

    # resize the image to look for interesting neighbours
    iter=0
    while True:
        startra,startdec=ra,dec
        tcopy=lt
        tcopy['dist']=np.sqrt((np.cos(dec*np.pi/180.0)*(tcopy['RA']-ra))**2.0+(tcopy['DEC']-dec)**2.0)*3600.0
        tcopy=tcopy[tcopy['dist']<180]
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
        if iter==10:
            break

    # now find the bounding box of the resulting collection
    ra,dec,size=find_bbox(tcopy)

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

    gals=Table.read(imagedir+'/tier1_ps1_wise_hetdex.fits')
    gals['raMean'].name='ra'
    gals['decMean'].name='dec'
    pg=gals[(np.abs(gals['ra']-ra)<size) & (np.abs(gals['dec']-dec)<size)]
    del(gals)

    gals=Table.read(imagedir+'/wise/allwise_HETDEX_full_radec.fits')
    pwg=gals[(np.abs(gals['ra']-ra)<size) & (np.abs(gals['dec']-dec)<size)]
    del(gals)

    ots=ot[(np.abs(ot['RA']-ra)<size) & (np.abs(ot['DEC']-dec)<size)]
    print ra,dec
    print ots['RA','DEC']

    pshdu=extract_subim(imagedir+'/downloads/'+psmaps[i],ra,dec,size,hduid=1)
    lhdu=extract_subim(lofarfile,ra,dec,size)
    firsthdu=extract_subim(imagedir+'/downloads/'+firstmaps[i],ra,dec,size)
    
    show_overlay(lhdu,pshdu,ra,dec,size,firsthdu=firsthdu,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,save_name=psimage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title)
    
    show_overlay(lhdu,pshdu,ra,dec,size,overlay_cat=ot,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,plotpos=[(pg,'x'),(pwg,'+')],show_lofar=False,save_name=pspimage,no_labels=True,title=title)
    whdu=extract_subim(imagedir+'/downloads/'+wisemaps[i],ra,dec,size)
    show_overlay(lhdu,whdu,ra,dec,size,firsthdu=firsthdu,overlay_cat=ot,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,save_name=wiseimage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,noisethresh=0)

    with open(manifestname,'w') as manifest:
        manifest.write('%i,%s,%s,%s,%s,%f,%f,%f\n' % (i,psimage,pspimage,wiseimage,sourcename,ra,dec,size*3600.0))

    os.system('mogrify -quality 90 -trim '+sourcename+'*.png')
