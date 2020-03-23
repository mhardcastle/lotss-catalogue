#!/usr/bin/env python3

'''
get_msource_flags
add flags to catalogue based on sub-flowchart of compact isolated m sources
'''

import os
import numpy as np

from astropy.table import Table, Column, join
import astropy.coordinates as ac
import astropy.units as u

from lofar_source_sorter_dr2 import Mask, Masks_disjoint_complete


#################################################################################



path = '/data2/wwilliams/projects/lofar_surveys/LoTSS-DR2-Feb2020/'
lofarcat_file_srt = path+'LoTSS_DR2_rolling.srl_0h.sorted_step1.fits'


lofarcat = Table.read(lofarcat_file_srt)



#################################################################################
# the following come from outputs from Lara/Philip for compact isolated m sources

#################################################################################
### msource1_flag
#0: no match
#1: accept ML of the source
#2: accept ML of the gaussian with highest ML
#3: deblend and accept both gaussians
#4: deblend workflow
#5: LOFAR galaxy zoo 

msource_cat_file = path+'msources/step1_isol_msources_flowchart_v1.fits'
msource_cat = Table.read(msource_cat_file)

if 'msource1_flag' in lofarcat.colnames:
    lofarcat.remove_column('msource1_flag')
if 'MC_flag1' in lofarcat.colnames:
    lofarcat.remove_column('MC_flag1')
lofarcat.sort('Source_Name')
tt=join(lofarcat, msource_cat, join_type='left', keys=['Source_Name'])
tt['M_Diagnosis_Code'].fill_value = -1
tt['MC_flag1'].fill_value = -1
tt = tt.filled()
tt.sort('Source_Name')
tt.rename_column('M_Diagnosis_Code','msource1_flag')
tt.rename_column('MC_flag1','MC1_flag1')


lofarcat.add_column(tt['msource1_flag'])
lofarcat.add_column(tt['MC1_flag1'])

# 
# compact non-isolated m sources
msource_cat_file = path+'msources/step1_nonisol_msources_flowchart_v1.fits'
msource_cat = Table.read(msource_cat_file)

if 'msource2_flag' in lofarcat.colnames:
    lofarcat.remove_column('msource2_flag')
if 'MC_flag2' in lofarcat.colnames:
    lofarcat.remove_column('MC_flag2')
lofarcat.sort('Source_Name')
tt=join(lofarcat, msource_cat, join_type='left', keys=['Source_Name'])
tt['M_Diagnosis_Code'].fill_value = -1
tt['MC_flag1'].fill_value = -1
tt = tt.filled()
tt.sort('Source_Name')
tt.rename_column('M_Diagnosis_Code','msource2_flag')
tt.rename_column('MC_flag1','MC2_flag1')


lofarcat.add_column(tt['msource2_flag'])
lofarcat.add_column(tt['MC2_flag1'])



#################################################################################

## write output file
lofarcat.write(lofarcat_file_srt, overwrite=True)
