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
from download_image_files import LofarMaps
import glob

scale=3600.0 # scaling factor for region sizes

def mogrify(filename):
    command='mogrify -quality 90 -trim '+filename
    p = subprocess.Popen(command,shell=True)
    return p


if __name__=='__main__':

    t=Table.read(glob.glob('sources_*.fits')[0])
    ot=Table.read(glob.glob('components_*.fits')[0])
    lm=LofarMaps(stay_in_imagedir=False)
    
    for i in range(0,len(t)):

        r=t[i]
        print(r)
        lofarfile=os.environ['IMAGEDIR']+'/'+lm.find(r['RA'],r['DEC'])
        sourcename=r['Source_Name']
        # find components of this source
        cfilter=ot['Parent_Source']==sourcename
        ctable=ot[cfilter]
        if len(ctable)==0:
            raise RuntimeError('Source has no components!')
        
        cimage=sourcename+'_C.png'
        cpimage=sourcename+'_Cp.png'

        if os.path.isfile(cimage):
            print('Selected output file exists already')
            continue

        ra,dec=r['RA'],r['DEC']

        marker_ra=None
        marker_dec=None
        title=None
        csize=r['LAS']*1.5
        if np.isnan(csize):
            csize=r['DC_Maj']*2*1.5
        
        # now find the bounding box
        ra,dec,size=find_bbox(ctable,scale=scale)
        
        if np.isnan(size) or size<csize:
            ra=r['RA']
            dec=r['DEC']
            size=csize

        print('size is',size)

        size/=3600.0

        ots=ctable

        ls=[]
        cs=[]
        for nr in ots:
            if nr['Parent_Source']==sourcename:
                ls.append('solid')
                cs.append('green')
            else:
                ls.append('dashed')
                cs.append('red')

        lhdu=extract_subim(lofarfile,ra,dec,size*2)
        lhdu.writeto(sourcename+'.fits')
        try:
            peak=r['Peak_flux']/1000
        except:
            peak=None

        if not np.isnan(r['optRA']):
            marker_ra=r['optRA']
            marker_dec=r['optDec']
        else:
            marker_ra=None
            marker_dec=None
            
        plist=[]
        error=False
        if not os.path.isfile(cimage):
            try:
                show_overlay(lhdu,lhdu,ra,dec,size,firsthdu=None,overlay_cat=ots,show_grid=False,overlay_scale=scale,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,save_name=cimage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,lw_ellipse=3,ellipse_style=ls,ellipse_color=cs,peak=peak,ppsize=350,cmap='inferno',lofar_colorscale=True,show_lofar=False)
                show_overlay(lhdu,lhdu,ra,dec,size,firsthdu=None,overlay_cat=None,show_grid=False,overlay_scale=scale,coords_color='red',coords_ra=None,coords_dec=None,coords_lw=3,lw=2,save_name=cpimage,no_labels=True,marker_ra=None,marker_dec=None,marker_lw=3,marker_color='cyan',title=title,lw_ellipse=3,ellipse_style=ls,ellipse_color=cs,peak=peak,ppsize=350,cmap='inferno',lofar_colorscale=True,show_lofar=False,plot_coords=False)
            except (Exception, RuntimeError) as e:
                error=True
            if error:
                print('*** image build failed! (%s) ***' % str(e))
                continue
            else:
                os.system('mogrify -quality 90 -trim '+sourcename+'*.png')
