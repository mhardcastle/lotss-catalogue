import string
import sys
import os
import numpy as np
from astropy.table import Table, Column, vstack




path = '/Users/w.williams/projects/lofar_surveys/DR2/'

ml_cat_file = path+'GradientBoostingClassifier_M2_31504_17F_TT42_B1_rd_dc_opt_tlv_high/pred_thresholds_full_13h.csv'
ml_cat = Table.read(ml_cat_file, format='csv')

h='13h'
#lofarcat_file_srt = path+'LoTSS_DR2_{version}.srl_{h}.lr-full.presort.fits'.format(version=version,h=h)
#scat = Table.read(lofarcat_file_srt)

# manually down for now ... but could be determined from np.unique
dup_names = ['ILTJ104642.23+371846.1',
'ILTJ113308.22+605332.4',
'ILTJ123019.45+564340.7',
'ILTJ125313.61+384929.8',
'ILTJ130927.32+500801.1',
'ILTJ163151.37+551815.7']

for dup in dup_names:
    si = np.where(ml_cat['Source_Name'] == dup)[0]
    
    i = 0
    for sii in si:
    
        ml_cat['Source_Name'][sii] = ml_cat['Source_Name'][sii]+string.ascii_letters[i]
        
        i += 1
    


#scat.write (path+'LoTSS_DR2_{version}.srl_{h}.lr-full.presort.fits'.format(version=version,h=h),overwrite=True)
ml_cat.write (path+'GradientBoostingClassifier_M2_31504_17F_TT42_B1_rd_dc_opt_tlv_high/pred_thresholds_full_13h_fixnames.fits',overwrite=True)
