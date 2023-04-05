#!/usr/bin/env python3

'''
make the selection of sources for Martin to feed into LGZ
weave priority , mostly flux above 4
'''
import sys
import os
import numpy as np

from astropy.table import Table, Column, join, vstack
import astropy.coordinates as ac
import astropy.units as u

#from lofar_source_sorter_dr2 import Mask, Masks_disjoint_complete


#################################################################################


mfix_field = 'P138+37'
#priority = '1'
#priority = '2'


version = 'v110'
path = '/Users/w.williams/projects/lofar_surveys/DR2/'
#lofarcat_file_srt = path+'LoTSS_DR2_{version}.srl_{h}.lr-full.sorted_step1.fits'.format(version=version,h=h)
lofarcat_file_srt = path+f'LoTSS_DR2_{version}.srl_13h.lr-full.sorted_step2_flux4.hdf5'


lofarcat = Table.read(lofarcat_file_srt)

for priority in ['mfix']:
    sel_file = path+f'lgz_selection_oct13_mfix/LoTSS_DR2_{version}.srl_13h.lr-full.sorted_step2_flux4.lgz_weave_selection_{priority}.fits'
    if priority in ['1' , '2' , '1a', '3']:
        sel_pri = (lofarcat['WEAVE_priority{i}'.format(i=priority)] == True)
    elif priority == '1b':
        sel_pri = (lofarcat['WEAVE_priority1'] == True) & (lofarcat['WEAVE_priority1a'] == False)  # special case of 1 and not 1a  make up 1b
    # fields not included in dr1 
    elif priority == 'mfix':
        sel_pri = (lofarcat['Mosaic_ID']==mfix_field)
    elif priority == 'hetdex':
        sel_pri = (lofarcat['hetdex_field']>0)
    # remaining non-priority areas
    elif priority == 'last':
        sel_pri = (~lofarcat['WEAVE_priority1']) & (~lofarcat['WEAVE_priority2']) & (~lofarcat['WEAVE_priority3']) & (lofarcat['hetdex_field']==0)

    sel =  sel_pri & (
            (lofarcat['FC_flag2'] == 5) | \
            (lofarcat['FC_flag2'] == 12) | \
            (lofarcat['FC_flag2'] == 15) | \
            (lofarcat['FC_flag2'] == 26) | \
            (lofarcat['FC_flag2'] == 20) )



    lofarcat1 = lofarcat[sel]

    lofarcat1.write(sel_file, overwrite=True)

    sel_file = path+f'lgz_selection_oct13_mfix/LoTSS_DR2_{version}.srl_13h.lr-full.sorted_step2_flux4.prefilter_lgz_weave_selection_{priority}.fits'

    sel = sel_pri & (
            (lofarcat['FC_flag2'] == 6) | \
            (lofarcat['FC_flag2'] == 13) | \
            (lofarcat['FC_flag2'] == 16) | \
            (lofarcat['FC_flag2'] == 27) | \
            (lofarcat['FC_flag2'] == 21) )



    lofarcat1 = lofarcat[sel]

    lofarcat1.write(sel_file, overwrite=True)
