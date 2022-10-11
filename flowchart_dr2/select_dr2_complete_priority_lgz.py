#!/usr/bin/env python3

'''
make the selection of sources for Martin to feed into LGZ
sections to compllete ids
'''
import sys
import os
import numpy as np

from astropy.table import Table, Column, join, vstack
import astropy.coordinates as ac
import astropy.units as u

from lofar_source_sorter_dr2 import Mask, Masks_disjoint_complete


#################################################################################


#priority = '1'
#priority = '2'


h = str(sys.argv[1])
if 'h' not in h:
    h+='h'
if h not in  ['0h','13h']:
    print('unknown field code (should be 0h or 13h)',h)
    sys.exit(1)
    
    
path = '/Users/w.williams/projects/lofar_surveys/DR2/'
lofarcat_file_srt = path+f'LoTSS_DR2_{version}.srl_{h}.lr-full.sorted_step2_flux4.hdf5'


lofarcat = Table.read(lofarcat_file_srt)

for priority in ['1']:
    sel_file_lgz = path+f'lgz_selection_apr13/LoTSS_DR2_{version}.srl_{h}.lr-full.sorted_step2_flux4.lgz_weave_selection_{priority}.fits'
    sel_file_pf = path+f'lgz_selection_apr13/LoTSS_DR2_{version}.srl_{h}.lr-full.sorted_step2_flux4.prefilter_lgz_weave_selection_{priority}.fits'
    
    sel_pri = (lofarcat['Complete_priority{i}'.format(i=priority)] == True)
        
    # the different fields have different flowcharts
    if h =='13h':
        FC_pri = [4, 9, 14, 17, 28]
    elif h =='0h':
        FC_pri = [4, 10, 17, 28]

    sel_FC_pri = np.zeros(len(lofarcat), dtype=bool)
    for pi in FC_pri:
        sel_FC_pri = sel_FC_pri | (lofarcat['FC_flag2'] == pi)
    
    # for LGZ select the TBD and ML=LGZ sources
    sel =  sel_pri &  (lofarcat['ML_flag'] == 0) & sel_FC_pri 
    lofarcat1 = lofarcat[sel]
    lofarcat1.write(sel_file_lgz, overwrite=True)

    # for prefilter select the TBD and ML=LR sources
    sel = sel_pri &  (lofarcat['ML_flag'] == 1) & sel_FC_pri
    lofarcat1 = lofarcat[sel]
    lofarcat1.write(sel_file_pf, overwrite=True)
