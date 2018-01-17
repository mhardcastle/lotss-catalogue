import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, join, Column
import astropy.units as u
import astropy.coordinates as ac
import utils.plot_util as pp
import os

from lofar_source_sorter import Mask


if __name__=='__main__':
    
    version = '0.8'

    path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'

    lofarcat_file_srt = path+'LOFAR_HBA_T1_DR1_catalog_v0.95_masked.srl.fixed.sorted.fits'

    # Source catalogue
    lofarcat = Table.read(lofarcat_file_srt)
    
    
    # LGZ output
    lgz_cat_file = os.path.join(path,'lgz_v1/HETDEX-LGZ-cat-v0.8-filtered-zooms.fits') 
    lgz_component_file = os.path.join(path,'lgz_v1/lgz_components.txt')


    lgz_components = Table.read(lgz_component_file, format='ascii', names=['lgz_component', 'lgz_src', 'lgz_flag'])
    
    
    
    
    
    
    selind = (lofarcat['FC_flag'] == 19)


    selindNN = (lofarcat['NN_idx'][selind])
    
    lofarcat_pairs = lofarcat[selind]
    lofarcat_pairsNN = lofarcat[selindNN]
    
    inlgz = []
    for name in lofarcat_pairs['Source_Name']:
        inlgz.append(name in lgz_components['lgz_component'])        
    inlgz = np.array(inlgz)
    
    print len(lofarcat_pairs)
    
    
    print 'excluding {} sources already associated in lgz'.format(np.sum(inlgz))
    lofarcat_pairs = lofarcat_pairs[~inlgz]
    lofarcat_pairsNN = lofarcat_pairsNN[~inlgz]
    
    
    inpair = np.zeros(len(lofarcat_pairs),dtype=bool)
    keep = np.ones(len(lofarcat_pairs),dtype=bool)
    for ni in range(len(lofarcat_pairs)):
        name = lofarcat_pairs['Source_Name'][ni]
        nameNN = lofarcat_pairsNN['Source_Name'][ni]
        
        if nameNN in lofarcat_pairs['Source_Name']:
            inpair[ni] = 1
            if keep[ni]:
                nni = np.where(lofarcat_pairs['Source_Name'] == nameNN)[0][0]
                keep[nni] = 0
        #print ni, name,  lofarcat_pairsNN['Source_Name'][ni], inpair[ni], keep[ni]
                 
       
    print np.sum(inpair)
    print np.sum(inpair&keep)
    print np.sum(~inpair)
    
    inlgz2 = (lofarcat_pairsNN['ID_flag']==3210)
    
    print 'excluding {} sources with NN already in lgz2'.format(np.sum(inlgz2))
    
    lofarcat_pairs = lofarcat_pairs[~inlgz2]
    lofarcat_pairsNN = lofarcat_pairsNN[~inlgz2]
    
    inpair = np.zeros(len(lofarcat_pairs),dtype=bool)
    keep = np.ones(len(lofarcat_pairs),dtype=bool)
    for ni in range(len(lofarcat_pairs)):
        name = lofarcat_pairs['Source_Name'][ni]
        nameNN = lofarcat_pairsNN['Source_Name'][ni]
        
        if nameNN in lofarcat_pairs['Source_Name']:
            inpair[ni] = 1
            if keep[ni]:
                nni = np.where(lofarcat_pairs['Source_Name'] == nameNN)[0][0]
                keep[nni] = 0
        #print ni, name,  lofarcat_pairsNN['Source_Name'][ni], inpair[ni], keep[ni]
                 
        
    print np.sum(inpair)
    print np.sum(inpair&keep)
    print np.sum(~inpair)
    
    selind = (inpair&keep) | (~inpair)
    selindNN = lofarcat_pairs[selind]['NN_idx']
    
    ra_pairs = 0.5*(lofarcat_pairs[selind]['RA'] + lofarcat[selindNN]['RA'])
    dec_pairs = 0.5*(lofarcat_pairs[selind]['DEC'] + lofarcat[selindNN]['DEC'])
    sep_pairs = lofarcat_pairs[selind]['NN_sep']
    
    lofarcat_pairs2 = lofarcat_pairs[selind]
    #lofarcat_pairs2['RA'] = ra_pairs
    #lofarcat_pairs2['DEC'] = dec_pairs
    
    lofarcat_pairs2.add_column(Column(2*sep_pairs,'Size'))
    lofarcat_pairs2.add_column(Column(ra_pairs,'Cent_RA'))
    lofarcat_pairs2.add_column(Column(dec_pairs,'Cent_DEC'))
    
    lofarcat_pairs2.write('doubles/sample_all_src_clean_small_nisol_nclustered_S_nlr_NNnlr_simflux_dist_lgz.fits', overwrite=True)
    
    
    lofarcat_pairs1 = lofarcat_pairs[~inpair]
    lofarcat_pairs1.write('doubles/sample_all_src_clean_small_nisol_nclustered_S_nlr_NNnlr_simflux_dist_lgz1.fits', overwrite=True)
