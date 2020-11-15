import numpy as np
import sys
import os
from astropy.table import Table
from astropy.io import fits
import glob
from subim import extract_subim
from overlay import find_noise_area,find_noise
from download_image_files import LofarMaps
from astropy.coordinates import SkyCoord
import astropy.units as u

# Script to generate cutouts and numpy arrays from LoTSS source catalogues for use with LoMorph (Mingo et al. 2019) and Barkus+ ridgeline host IDs code. Takes a LoTSS source catalogue, and writes out a directory of fits cutouts, and npy array thresholded at 4 x rms with dynamic range threshold applied, as described in Mingo+19. Written by J.Croston.

def get_fits(fra,fdec,fsource,fsize):

    sc=SkyCoord(fra*u.deg,fdec*u.deg,frame='icrs')
    s=sc.to_string(style='hmsdms',sep='',precision=2)
    name=fsource
    newsize=2.5*fsize/3600.0

    lm=LofarMaps()
    print fra,fdec
    mosname=lm.find(fra,fdec)
    filename=os.environ['IMAGEDIR']+'/'+mosname
    hdu=extract_subim(filename,fra,fdec,newsize)
    if hdu is not None:
        hdu.writeto('cutouts/'+name+'.fits',overwrite=True)
        flag=0
    else:
        print 'Cutout failed for '+str(fsource)
        flag=1

    return flag

# Read in tables

sourcecat = sys.argv[1]
summary=Table.read(sourcecat)

# Go through summary table, use LGZ_Size if present, else Maj to define cutout reg. 

for row in summary:
    maj=row['Maj']
    lgzsiz=row['LGZ_Size']
    ssource=row['Source_Name']
    #ssize=row['Size']
    ssource=ssource.rstrip()
    print "source is "+ssource
    sra=row['RA']
    sdec=row['DEC']
    flux=row['Peak_flux']
    rms=row['Isl_rms']
    
    if np.isnan(maj):
        ssize=lgzsiz
    else:
        ssize=maj
   
    # Take 20 arcsec as minimum cutout size
    if ssize<20.0:
        ssize=20.0

    flag=get_fits(sra,sdec,ssource,ssize)
    # make output locations if needed
    if os.path.isdir(cutouts):
        continue
    else:
        os.mkdir('cutouts')
    if os.path.isdir(4rms):
        continue
    else:
        os.mkdir('4rms')
    
    cutout='cutouts/'+ssource+'.fits'
    lmd=LofarMaps()
    dname=os.environ['IMAGEDIR']+'/'+lmd.find(sra,sdec)
    if os.path.isfile(dname):
        lhdu=fits.open(dname)
        if flag==0:
            nlhdu=fits.open(cutout)

            # Write out npy arrays applying 4xrms and dyn range thresholds (cf Mingo+19)
            
            d=nlhdu[0].data
            peak=flux
            dyncut=50.0
            ratio=peak/dyncut
            print peak,rms,ratio
            if rms<ratio:
                print "rms < ratio"
                thres=ratio

            else:
                print "rms > ratio"
                print rms
                thres=4.0*rms
 
            d[d<thres]=np.nan
            np.save("4rms/"+ssource+".npy",d)

