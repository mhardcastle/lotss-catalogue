#!/usr/bin/env python3

'''
get_msource_flags
add flags to catalogue based on sub-flowchart of compact isolated m sources
'''
import sys
import os
import numpy as np

from astropy.table import Table, Column, join
import astropy.coordinates as ac
import astropy.units as u

from lofar_source_sorter_dr2 import Mask, Masks_disjoint_complete


#################################################################################





if len(sys.argv) == 1:
    print("Usage is : python get_ml_flags.py field_code ")
    print('E.g.: python get_ml_flags.py 0 ')
    sys.exit(1)

h = str(sys.argv[1])
if 'h' not in h:
    h+='h'
if h not in  ['0h','13h','n0h','n13h','s0h','s13h']:
    print('unknown field code (should be 0h or 13h)',h)
    sys.exit(1)

path = '/data2/wwilliams/projects/lofar_surveys/LoTSS-DR2-Feb2020/'
lofarcat_file_srt = path+'LoTSS_DR2_v100.srl_{h}.lr-full.sorted_step1.fits'.format(h=h)



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
