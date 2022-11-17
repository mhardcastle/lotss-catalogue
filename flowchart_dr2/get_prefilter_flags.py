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
if h not in  ['0h','13h']:
    print('unknown field code (should be 0h or 13h)',h)
    sys.exit(1)

path = '/Users/w.williams/projects/lofar_surveys/DR2/'
lofarcat_file_srt = path+'LoTSS_DR2_{version}.srl_{h}.lr-full.sorted_step2_flux4.hdf5'.format(version=version,h=h)
#lofarcat_file_srt = path+'LoTSS_DR2_{version}.srl_{h}.lr-full.presort.hdf5'.format(version=version,h=h)



#################################################################################
# the following come from outputs from flow_python/extract_prefilter_out.py (run on herts)
p_cat_file = path+'LoTSS_DR2_{version}.srl_{h}.prefilter_outputs.fits'.format(version=version,h=h)
p_cat = Table.read(p_cat_file)

lofarcat = Table.read(lofarcat_file_srt)
lofarcat.sort('Source_Name')

if 'Prefilter' in lofarcat.colnames:
    lofarcat.remove_columns('Prefilter')

for i in range(len(lofarcat)):
    lofarcat['Source_Name'][i] = lofarcat['Source_Name'][i].strip()

tt = join(lofarcat, p_cat, keys='Source_Name')
tt.sort('Source_Name')

lofarcat.add_column(tt['Prefilter'])


## write output file
lofarcat.write(lofarcat_file_srt, overwrite=True, serialize_meta=True)
