#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
from astropy.table import Table,vstack
import numpy as np
import sys
import os
import glob
from subim import extract_subim
from overlay import show_overlay

scale=3600.0 # scaling factor for region sizes

def get_mosaic_name(name):
    globst=os.environ['IMAGEDIR']+'/mosaics/'+name.rstrip()+'*'
    print 'Looking for',globst
    g=glob.glob(globst)
    if len(g)>0:
        return g[0]
    else:
        raise RuntimeError('No mosaic called '+name)

if __name__=='__main__':

    imagedir=os.environ['IMAGEDIR']
    
    # all LOFAR sources
    ot=Table.read(imagedir+'/LOFAR_HBA_T1_DR1_catalog_v0.2.srl.fits')

    # read lists
    lines=[l.rstrip().split() for l in open('maplist.txt').readlines()]
    lofarmaps={}
    psmaps={}
    firstmaps={}
    for l in lines:
        name=l[0]
        lofarmaps[name]=l[1]
        psmaps[name]=l[2]
        firstmaps[name]=l[4]
    
    # read association check list

    check_sources=[]
    check_positions={}
    lines=open('assoc-check.txt').readlines()
    for l in lines:
        bits=l.rstrip().split()
        name=bits[0]
        if name not in check_sources:
            check_sources.append(name)
            check_positions[name]=[]
        check_positions[name].append((float(bits[1]),float(bits[2])))
    
    # read manifest for positions and sizes
    lines=open('combined_manifest.csv').readlines()
    image={}
    for l in lines[1:]:
        bits=l.split(',')
        name=bits[4]
        image[name]=(float(bits[5]),float(bits[6]),float(bits[7]))
    
    # read external tables

    #gals=Table.read(imagedir+'/tier1_ps1_wise_hetdex.fits')
    #gals['raMean'].name='ra'
    #gals['decMean'].name='dec'
    #wgals=Table.read(imagedir+'/wise/allwise_HETDEX_full_radec.fits')

    # now go through the sources
    for name in check_sources:
        print lofarmaps[name],psmaps[name],firstmaps[name]
        outname=name+'.png'
        if os.path.isfile(outname):
            print name,'output file exists, skipping'
            continue

        print '*** Doing source',name,'***'

        lofarfile=os.environ['IMAGEDIR']+'/'+lofarmaps[name]

        ra,dec,size=image[name]
        size/=3600.0

        #pg=gals[(np.abs(gals['ra']-ra)<size) & (np.abs(gals['dec']-dec)<size)]

        #pwg=wgals[(np.abs(wgals['ra']-ra)<size) & (np.abs(wgals['dec']-dec)<size)]

        pshdu=extract_subim(imagedir+'/downloads/'+psmaps[name],ra,dec,size,hduid=1)
        lhdu=extract_subim(lofarfile,ra,dec,size)
        firsthdu=extract_subim(imagedir+'/downloads/'+firstmaps[name],ra,dec,size)
        plot_positions=check_positions[name]
        print 'Plotting positions',plot_positions
        plot_ra=[]
        plot_dec=[]
        for p in plot_positions:
            plot_ra.append(p[0])
            plot_dec.append(p[1])
        show_overlay(lhdu,pshdu,ra,dec,size,firsthdu=firsthdu,overlay_cat=ot,overlay_scale=scale,lw=2,save_name=outname,no_labels=True,coords_ra=plot_ra,coords_dec=plot_dec)

