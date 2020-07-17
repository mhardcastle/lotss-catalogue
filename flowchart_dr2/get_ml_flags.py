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
# the following come from outputs from Lara/Philip for ML classification

#################################################################################

ml_cat_file = path+'GradientBoostingClassifier_lotss_dr2_0h_pred_thresholds.fits'
ml_cat = Table.read(ml_cat_file)

if 'msource1_flag' in lofarcat.colnames:
    lofarcat.remove_column('msource1_flag')
if 'MC_flag1' in lofarcat.colnames:
    lofarcat.remove_column('MC_flag1')
lofarcat.sort('Source_Name')
tt=join(lofarcat, ml_cat, join_type='left', keys=['Source_Name'])
tt['0.12'].fill_value = -1
tt = tt.filled()
tt.sort('Source_Name')
tt.rename_column('0.12','ML_flag')  # 1 for LR, 0 for LGZ


lofarcat.add_column(tt['ML_flag'])


#sys.exit()
#################################################################################

## write output file
lofarcat.write(lofarcat_file_srt, overwrite=True)
