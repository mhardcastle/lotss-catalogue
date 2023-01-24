#!/usr/bin/env python3

'''
make the selection of sources for Martin to feed into LGZ
the Msources with ML=LGZ that I thought had been done in LGZ but weren't
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
selno = '3'
version = 'v110'

path = '/Users/w.williams/projects/lofar_surveys/DR2/'
#lofarcat_file_srt = path+'LoTSS_DR2_{version}.srl_{h}.lr-full.sorted_step1.fits'.format(version=version,h=h)
lofarcat_file_srt = path+'LoTSS_DR2_{version}.srl_0h.lr-full.sorted_step3_flux4.hdf5'.format(version=version)
sel_file = path+'lgz_selection_2023jan19/LoTSS_DR2_{version}.srl_0h.lr-full.sorted_step2_flux4.lgz_selection_{i}.fits'.format(i=selno,version=version)


lofarcat = Table.read(lofarcat_file_srt)


sel_pri = (lofarcat['WEAVE_priority{i}'.format(i=priority)] == True)

sel =  sel_pri & (lofarcat['ML_flag']==0) & (
        (lofarcat['FC_flag3'] == 19) | \
        (lofarcat['FC_flag3'] == 29)| \
        (lofarcat['FC_flag3'] == 41) )

lofarcat1 = lofarcat[sel]

lofarcat1.write(sel_file, overwrite=True)

sel_file = path+'lgz_selection_2023jan19/LoTSS_DR2_{version}.srl_0h.lr-full.sorted_step2_flux4.prefilter_lgz_selection_{i}.fits'.format(i=selno,version=version)

# note that the ML_flag==0 selection is already made in the flowchart by selecting those without prefilter outcomes
sel =  sel_pri & (
        (lofarcat['FC_flag3'] == 27) | \
        (lofarcat['FC_flag3'] == 37) | \
        (lofarcat['FC_flag3'] == 49))




lofarcat1 = lofarcat[sel]

lofarcat1.write(sel_file, overwrite=True)



sel_file = path+'lgz_selection_2023jan19/LoTSS_DR2_{version}.srl_0h.lr-full.sorted_step2_flux4.deblend_lgz_selection_{i}.fits'.format(i=selno,version=version)

# note have ML_flag==0 selection been done  -- TBC 
sel =  sel_pri & (
        (lofarcat['FC_flag3'] == 18) | \
        (lofarcat['FC_flag3'] == 40))




lofarcat1 = lofarcat[sel]

lofarcat1.write(sel_file, overwrite=True)
