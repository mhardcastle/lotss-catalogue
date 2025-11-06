# Make number counts

from astropy.table import Table
import numpy as np
from tqdm import tqdm
import healpy as hp
import matplotlib.pyplot as plt
import os
import pickle
import sys

NSIDE=16 # pixels to work with

# Read the area map

amap=hp.fitsfunc.read_map('area_hp.fits')
bmap=np.where(amap<1180,amap,1181)
bmap[np.isnan(amap)]=0
blow=hp.nside2pixarea(NSIDE,degrees=True)*hp.pixelfunc.ud_grade(bmap,NSIDE)/1181 # now in deg^2

# Read the rms map

rmap=hp.fitsfunc.read_map('rms_hp.fits')
rmap[np.isnan(rmap)]=hp.UNSEEN
rlow=hp.pixelfunc.ud_grade(rmap,NSIDE)

# Read the source cat

ti=Table.read('LoTSS_DR3_v0.4.srl.fits')
ti['pixel']=hp.ang2pix(NSIDE,ti['RA'],ti['DEC'],lonlat=True)
ti['Total_flux']/=1000 # Jy
up=sorted(np.unique(ti['pixel']))

for pix in tqdm(up):
    t=ti[ti['pixel']==pix]
    picklefile=f'sourcecount-{pix}.npy'
    if os.path.isfile(picklefile):
        print('loading pickle')
        with open(picklefile,'rb') as infile:
            bc,bins,dhist,hist=pickle.load(infile)
            bs=bins[1:]-bins[:-1]
            bs=bs[:len(dhist)]
    else:
        area=blow[pix]
        rms=rlow[pix]
        if rms<0 or area<7: continue
        print('Pixel',pix)
        print('Area',area,'deg^2')
        print('rms',rms,'Jy/beam')
        cutoff=rms*50
        print('Cutoff will be',cutoff,'Jy')
        #cutoff=1e-3
        maxflux=np.max(t['Total_flux'])
        t=t[t['Total_flux']>cutoff]
        print('Number of sources after completeness cut is',len(t))

        bins=np.logspace(np.log10(cutoff),np.log10(maxflux)*1.01,20)
        cbins=0.5*(bins[:-1]+bins[1:])
        ds=bins[1:]-bins[:-1]
        bc=(bins[1:]+bins[:-1])/2
        bs=bins[1:]-bins[:-1]
        hist,_=np.histogram(t['Total_flux'],bins=bins)
        dhist=hist/area
        cutoff=np.where(hist==0)[0]
        if len(cutoff)>0:
            print('cutting at',cutoff[0])
            bc=bc[:cutoff[0]]
            dhist=dhist[:cutoff[0]]
            bs=bs[:cutoff[0]]
        with open(picklefile,'wb') as outfile:
            pickle.dump([bc,bins,dhist,hist],outfile)
    plt.plot(bc,dhist/bs,color='red',alpha=0.03,lw=1)

if len(sys.argv)>1:
    with open('median.npy','rb') as infile:
        median_bc,median=pickle.load(infile)
    plt.plot(median_bc,median,color='blue',lw=3)
    
plt.xlabel('Total flux (Jy)')
plt.ylabel('$N(S)$ Jy$^{-1}$ deg$^{-2}$')
plt.xscale('log')
plt.yscale('log')

if len(sys.argv)>2:
    plt.savefig('sourcecounts.pdf')
else:
    plt.show()
