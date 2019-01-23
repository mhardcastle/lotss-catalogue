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

scale=3600.0 # scaling factor for region sizes

if __name__=='__main__':

    imagedir=os.environ['IMAGEDIR']
    tname=sys.argv[1]
    lname=sys.argv[1].replace('.fits','-list.txt')
    t=Table.read(tname)
    # Annotated PyBDSF table
    ot=Table.read(os.environ['LOTSS_COMPONENT_CATALOGUE'])
    # large source table for neighbours
    lt=ot[(ot['Total_flux']>3) & (ot['Maj']>8)]

    # read lists
    lines=[l.rstrip().split() for l in open(lname).readlines()]
    lofarmaps=[l[1] for l in lines]
    psmaps=[l[2] for l in lines]
    firstmaps=[l[4] for l in lines]

    start=int(sys.argv[2])
    try:
        end=int(sys.argv[3])+1
    except:
        end=start+1
    for i in range(start,end):
        r=t[i]
        try:
            sourcename=r['Source_Name']
        except KeyError:
            # try to fix up old-format files
            t['Source_id'].name='Source_Name'
            sourcename=r['Source_Name']
        sourcename=sourcename.rstrip()
        print i,sourcename

        psimage=sourcename+'_PS.png'
        pspimage=sourcename+'_PSp.png'
        manifestname=sourcename+'-manifest.txt'
        print lofarmaps[i],psmaps[i]

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
        if '.fits' not in lofarfile:
            raise RuntimeError('Could not find mosaic fits file')
        
        if os.path.isfile(manifestname):
            print 'Selected output file exists already'
            continue

        from overlay import show_overlay

        ra,dec=r['RA'],r['DEC']

        try:
            marker_ra=r['ra']
            marker_dec=r['dec']
        except:
            marker_ra=None
            marker_dec=None

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
            size=300.0
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

        ots=ots[ots['Source_Name']!=""] # removes artefacts
        ls=[]
        for nr in ots:
            if nr['Component_Name']==r['Source_Name']:
                ls.append('solid')
            else:
                ls.append('dashed')

        pshdu=fits.open(imagedir+'/downloads/'+psmaps[i])
        lhdu=extract_subim(lofarfile,ra,dec,size)
        lhdu.writeto('lofar.fits',clobber=True)
        firsthdu=extract_subim(imagedir+'/downloads/'+firstmaps[i],ra,dec,size)
        try:
            peak==r['Peak_flux']/1000.0
        except:
            peak=None


        show_overlay(lhdu,pshdu,ra,dec,size,firsthdu=None,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_lw=3,lw=1,save_name=psimage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,peak=peak,plot_coords=False,show_grid=False,lw_ellipse=3,ellipse_style=ls,noisethresh=1.5)
        
        show_overlay(lhdu,pshdu,ra,dec,size,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_lw=3,lw=2,show_lofar=False,save_name=pspimage,no_labels=True,title=title,peak=peak,plot_coords=False,show_grid=False,lw_ellipse=3,ellipse_style=ls,noisethresh=1.5)

        with open(manifestname,'w') as manifest:
            manifest.write('%i,%s,%s,%s,%f,%f,%f\n' % (i,psimage,pspimage,sourcename,ra,dec,size*3600.0))

        os.system('mogrify -quality 90 -trim '+sourcename+'*.png')
