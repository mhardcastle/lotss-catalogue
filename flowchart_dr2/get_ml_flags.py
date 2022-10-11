#!/usr/bin/env python3

'''
get_msource_flags
add flags to catalogue based on sub-flowchart of compact isolated m sources
'''
import sys
import os
import numpy as np

from astropy.table import Table, Column, join, vstack
import astropy.coordinates as ac
import astropy.units as u

#from lofar_source_sorter_dr2 import Mask, Masks_disjoint_complete


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

path = '/Users/w.williams/projects/lofar_surveys/DR2/'
#lofarcat_file_srt = path+'LoTSS_DR2_{version}.srl_{h}.lr-full.sorted_step1.fits'.format(version=version,h=h)
lofarcat_file_srt = path+'LoTSS_DR2_{version}.srl_{h}.lr-full.presort.hdf5'.format(version=version,h=h)



#################################################################################
# the following come from outputs from Lara/Philip for ML classification

#################################################################################

if h == '0h':
    ml_cat_file = path+'GradientBoostingClassifier_lotss_dr2_0h_pred_thresholds.fits'
    ml_cat = Table.read(ml_cat_file)
    
    lofarcat = Table.read(lofarcat_file_srt)
    lofarcat.sort('Source_Name')

    ml_cat.rename_column('0.12','ML_flag')  # 1 for LR, 0 for LGZ
    tt = join(lofarcat, ml_cat)
    tt.sort('Source_Name')
    
    if 'ML_flag' in lofarcat.colnames:
        lofarcat['ML_flag'] = tt['ML_flag']
    else:
        lofarcat.add_column(tt['ML_flag'])


    ## write output file
    lofarcat.write(lofarcat_file_srt, overwrite=True, serialize_meta=True)

elif h == '13h':
    
    # concat in topcat and keep only 0.12 col
    ml_cat_file = path+'GradientBoostingClassifier_M2_31504_17F_TT42_B1_rd_dc_opt_tlv_high/pred_thresholds_full_13h_fixnames.fits'
    ml_cat = Table.read(ml_cat_file)

    lofarcat = Table.read(lofarcat_file_srt)
    lofarcat.sort('Source_Name')
    for i in range(len(lofarcat)):
        lofarcat['Source_Name'][i] = lofarcat['Source_Name'][i].strip()
    for i in range(len(ml_cat)):
        ml_cat['Source_Name'][i] = ml_cat['Source_Name'][i].strip()

    ml_cat.rename_column('0.12','ML_flag')  # 1 for LR, 0 for LGZ
    tt = join(lofarcat, ml_cat, keys='Source_Name')
    tt.sort('Source_Name')
    
    if 'ML_flag' in lofarcat.colnames:
        lofarcat['ML_flag'] = tt['ML_flag']
    else:
        lofarcat.add_column(tt['ML_flag'])


    ## write output file
    lofarcat.write(lofarcat_file_srt, overwrite=True, serialize_meta=True)

    #cmd = 'stilts tmatch2  in1={fin1} in2={fin2} ifmt2=csv matcher=exact values1=Source_Name values2=Source_Name join=all1 out=match.fits'.format(fin1=lofarcat_file_srt, fin2=ml_cat_file)
    #os.system(cmd)

    #tt = Table.read('match.fits',memmap=True)
    ##tt.keep_columns(['Source_Name_1','0.12'])
    #if 'ML_flag' in tt.colnames:
        #tt.remove_column('ML_flag')
    #tt.rename_column('0.12','ML_flag')  # 1 for LR, 0 for LGZ
    #tt.rename_column('Source_Name_1','Source_Name')  # 1 for LR, 0 for LGZ
    #tt.remove_column('Source_Name_2')


    ### write output file
    #tt.write(lofarcat_file_srt, overwrite=True, serialize_meta=True)

#sys.exit()
#################################################################################

