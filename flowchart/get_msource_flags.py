#!/usr/bin/python

'''
get_msource_flags
add flags to catalogue based on sub-flowchart of compact isolated m sources
'''

from lofar_source_sorter import Mask, Masks_disjoint_complete
import numpy as np
from astropy.table import Table, Column, join
import astropy.coordinates as ac
import astropy.units as u
import os

#################################################################################

path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'
lofarcat_file_srt = path+'LOFAR_HBA_T1_DR1_catalog_v0.95_masked.srl.fixed.presort.fits'



lofarcat = Table.read(lofarcat_file_srt)



#################################################################################
# the following come from outputs from Lara/Philip for compact isolated m sources

#################################################################################
### nhuge_2masx
#0: no match
#1: accept ML of the source
#2: accept ML of the gaussian with highest ML
#3: deblend and accept both gaussians
#4: deblend workflow
#5: LOFAR galaxy zoo 

msource_cat_file = path+'msources/LOFAR_flowchart.fixed.fits'
msource_cat = Table.read(msource_cat_file)

if 'msource1_flag' in lofarcat.colnames:
    lofarcat.remove_column('msource1_flag')
lofarcat.sort('Source_Name')
tt=join(lofarcat, msource_cat, join_type='left')
tt['Diagnosis_Code'].fill_value = -1
tt = tt.filled()
tt.sort('Source_Name')
tt.rename_column('Diagnosis_Code','msource1_flag')


lofarcat.add_column(tt['msource1_flag'])

#################################################################################

## write output file
lofarcat.write(lofarcat_file_srt, overwrite=True)