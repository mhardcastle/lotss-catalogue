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

from lofar_source_sorter_dr2 import Mask, Masks_disjoint_complete


#################################################################################


priority = '1'
#priority = '2'


path = '/data2/wwilliams/projects/lofar_surveys/LoTSS-DR2-Feb2020/'
#lofarcat_file_srt = path+'LoTSS_DR2_v100.srl_{h}.lr-full.sorted_step1.fits'.format(h=h)
lofarcat_file_srt = path+'LoTSS_DR2_v100.srl_13h.lr-full.sorted_step2_flux4.hdf5'


lofarcat = Table.read(lofarcat_file_srt)

for priority in ['1','2']:
    sel_file = path+'lgz_selection_nov18/LoTSS_DR2_v100.srl_13h.lr-full.sorted_step2_flux4.lgz_weave_selection_{i}.fits'.format(i=priority)
    if priority == '1' or priority == '2' or priority == '1a':
        sel_pri = (lofarcat['WEAVE_priority{i}'.format(i=priority)] == True)
    elif priority == '1b':
        sel_pri = (lofarcat['WEAVE_priority1'] == True) & (lofarcat['WEAVE_priority1a'] == False)  # special case of 1 and not 1a  make up 1b

    sel =  sel_pri & (
            (lofarcat['FC_flag2'] == 5) | \
            (lofarcat['FC_flag2'] == 12) | \
            (lofarcat['FC_flag2'] == 15) | \
            (lofarcat['FC_flag2'] == 26) | \
            (lofarcat['FC_flag2'] == 20) )



    lofarcat1 = lofarcat[sel]

    lofarcat1.write(sel_file, overwrite=True)

    sel_file = path+'lgz_selection_nov18/LoTSS_DR2_v100.srl_13h.lr-full.sorted_step2_flux4.prefilter_lgz_weave_selection_{i}.fits'.format(i=priority)

    sel = sel_pri & (
            (lofarcat['FC_flag2'] == 6) | \
            (lofarcat['FC_flag2'] == 13) | \
            (lofarcat['FC_flag2'] == 16) | \
            (lofarcat['FC_flag2'] == 27) | \
            (lofarcat['FC_flag2'] == 21) )



    lofarcat1 = lofarcat[sel]

    lofarcat1.write(sel_file, overwrite=True)
