#!/usr/bin/python

'''
flag_edge_sources.py
check all sources for proximity to edge (by presence of nan pixels within Maj), further check these sources for actual overlap with edge using presence of nanpixels within mask defined by the pybdsf ellipse shape. 
'''

import matplotlib
matplotlib.use('Agg')
from astropy.table import Table,vstack,Column
from astropy.io import fits
from astropy.wcs import WCS
import astropy.coordinates as ac
import astropy.units as u
import numpy as np
import sys
import os
import glob
from subim import extract_subim
import multiprocessing as mp 
import itertools
import time
import matplotlib.pyplot as plt
import pyregion

IMAGEDIR='/data/lofar/mjh/hetdex_v4'


def get_mosaic_name(name):
    if name == 'P21':
        name = 'P21-mosaic'
    g=glob.glob(IMAGEDIR+'/mosaics/'+name+'*')
    if len(g)>0:
        return g[0]
    else:
        raise RuntimeError('No mosaic called '+name)

    
    
def check_flagged(mosaicname,ra,dec,size):
    
    mosaicfile = get_mosaic_name(mosaicname)
    
    size = np.max((size, 0.001))
    lhdu = extract_subim(mosaicfile,ra,dec,size,verbose=False)
    
    if lhdu is None:
        print mosaicname,ra,dec,size
        return True
    
    data = lhdu[0].data
    
    flagged = np.any(np.isnan(data))
        
    return flagged


def check_flagged_region(mosaicname,ra,dec,maj,min,pa,plot=False):
    
    mosaicfile = get_mosaic_name(mosaicname)
    
    maj = np.max((maj, 0.001))
    lhdu = extract_subim(mosaicfile,ra,dec,2*maj,verbose=False)
    
    if lhdu is None:
        print mosaicname,ra,dec,maj
        return True
        
    region = pyregion.parse('fk5;ellipse({ra},{dec},{a},{b},{pa})'.format(ra=ra,dec=dec,a=maj,b=min,pa=pa-90))
    
    data = lhdu[0].data
    region_mask = region.get_mask(lhdu[0])
    
    flagged = np.any(np.isnan(data[region_mask]))
    
    if plot:
        plt.figure()
        plt.subplot(131)
        plt.imshow(data)
        plt.subplot(132)
        plt.imshow(region_mask)
        plt.subplot(133)
        plt.imshow(data*region_mask)
    
    return flagged


def check_flagged_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return check_flagged(*a_b)
    
if __name__ == '__main__':
    lofarcat_file = '/data/lofar/wwilliams/hetdex/LoTSS-DR1-July21-2017/LOFAR_HBA_T1_DR1_catalog_v0.95_masked.srl.fits'
    lofarcat_flagged_file = '/data/lofar/wwilliams/hetdex/LoTSS-DR1-July21-2017/LOFAR_HBA_T1_DR1_catalog_v0.95_masked.srl.edgeflags.fits'
    
    lofarcat_all = Table.read(lofarcat_file)
    
    ncpu = 16
    pool = mp.Pool(processes=ncpu)
    
    #for some reason this stops midway?? do it in chunks
    
    Ns = ncpu
    N = ncpu*1200
    nanflags = np.zeros(len(lofarcat_all),dtype=bool)
    for i in range(Ns+1):
        i2 = np.min((len(lofarcat_all),(i+1)*N)) 
        print 'chunk {}: {}-{} '.format(i,i*N, i2)
        starttime = time.time()
        lofarcat = lofarcat_all[i*N:i2]
        results = pool.map(check_flagged_star, itertools.izip(lofarcat['Mosaic_ID'],lofarcat['RA'],lofarcat['DEC'],lofarcat['Maj']/3600.))
        nanflags[i*N:i2] = np.array(results)
        endtime = time.time()
    
        print 'Took {:.2f} seconds for {} sources' .format(endtime-starttime, len(lofarcat))
    
    # take a closer look at the ones that have been flagged as having a nanpixel within Maj
    # apply region mask and check for nan pixels inside that only
    nanflags2 = nanflags.copy()
    
    for ti in np.where(nanflags)[0]:
        tt = lofarcat_all[ti]
        newflag = check_flagged_region(tt['Mosaic_ID'],tt['RA'],tt['DEC'],tt['Maj']/3600.,tt['Min']/3600.,tt['PA'])
        print ti, tt['Source_Name'], newflag
        nanflags2[ti] = newflag
    
    # finally check if any flagged sources overlap with other unflagged sources - these are probably bad too, unless it is a huge source that is flagged
    # could be more sophisticated and check for actual overlap, not merely proximity based on Maj... but this seems to work.
    nanflags3 = nanflags2.copy()
    
    
    catc = ac.SkyCoord(lofarcat_all['RA'],lofarcat_all['DEC'],unit='deg')
    for ti in np.where(nanflags2)[0]:
        tt = lofarcat_all[ti]
        if tt['Maj'] < 45:
            catt = ac.SkyCoord(tt['RA'],tt['DEC'],unit='deg')
            sep = catt.separation(catc)
            nearby = np.where((sep<tt['Maj']*u.arcsec) & (sep>0*u.arcsec))[0]
            if len(nearby) >=  1:
                for ni in nearby:
                    if not nanflags2[ni]:
                        print tt['Source_Name'], tt['Maj'], '- also flagging ', lofarcat_all['Source_Name'][ni]
                    
                        nanflags3[ni] = True
            #lofarcat_all[nearby]

    
    
    lofarcat_all.add_column(Column(nanflags, 'Edge_flag1'))
    lofarcat_all.add_column(Column(nanflags2, 'Edge_flag2'))
    lofarcat_all.add_column(Column(nanflags3, 'Edge_flag3'))
    
    lofarcat_all.keep_columns(['Source_Name', 'Edge_flag1','Edge_flag2','Edge_flag3'])
    lofarcat_all.write(lofarcat_flagged_file,overwrite=True)
    
    print '{} edge1 flags'.format(np.sum(lofarcat_all['Edge_flag1']))
    print '{} edge2 flags'.format(np.sum(lofarcat_all['Edge_flag2']))
    print '{} edge3 flags'.format(np.sum(lofarcat_all['Edge_flag3']))
    
    
    # plot where the edges are
    plt.figure()
    plt.scatter(lofarcat_all['RA'],lofarcat_all['DEC'],edgecolor='none')
    plt.scatter(lofarcat_all['RA'][nanflags],lofarcat_all['DEC'][nanflags],edgecolor='none')
    plt.scatter(lofarcat_all['RA'],lofarcat_all['DEC'],edgecolor='none')
    plt.scatter(lofarcat_all['RA'][nanflags2],lofarcat_all['DEC'][nanflags2],edgecolor='none')
    plt.scatter(lofarcat_all['RA'],lofarcat_all['DEC'],edgecolor='none')
    plt.scatter(lofarcat_all['RA'][nanflags3],lofarcat_all['DEC'][nanflags3],edgecolor='none')
    #plt.scatter(lofarcat['RA'][nanflags],lofarcat['DEC'][nanflags],c='y',edgecolor='none')
    plt.savefig('edgeflags.png')
    
    
    # check the ones highlighted in the issue
    checklist = ['ILTJ151524.99+543046.0',
                    'ILTJ150319.52+454944.5',
                    'ILTJ142450.82+502655.3',
                    'ILTJ142342.77+524831.0',
                    'ILTJ141236.02+512257.8',
                    'ILTJ124236.50+562844.6',
                    'ILTJ114646.72+562353.2',
                    'ILTJ122510.05+460000.4',
                    'ILTJ122510.79+455955.5',
                    'ILTJ125725.79+561852.2',
                    'ILTJ125726.89+561842.4']
    
    for n in checklist: 
        print n, n in lofarcat_all[nanflags]['Source_Name'] , n in lofarcat_all[nanflags2]['Source_Name'] , n in lofarcat_all[nanflags3]['Source_Name']
        
    ## hmm, two not flagged - one is in a messy part of the sky, the other is near the edge but doesn't quite make the cut


    ## testing
    
    nn = (lofarcat_all['Source_Name'] == 'ILTJ122510.05+460000.4') | (lofarcat_all['Source_Name'] == 'ILTJ122510.79+455955.5') 
    for ti in np.where(nn)[0]:
        tt = lofarcat_all[ti]
        flag = check_flagged(tt['Mosaic_ID'],tt['RA'],tt['DEC'],tt['Maj']/3600.)
        newflag =  check_flagged_region(tt['Mosaic_ID'],tt['RA'],tt['DEC'],tt['Maj']/3600.,tt['Min']/3600.,tt['PA'], plot=True)
        print ti, tt['Source_Name'], flag, newflag
