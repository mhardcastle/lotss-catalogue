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

from lofar_source_sorter_dr2 import Mask, Masks_disjoint_complete


#################################################################################





if len(sys.argv) == 1:
    print("Usage is : python get_artefact_flags.py field_code ")
    print('E.g.: python get_artefact_flags.py 0 ')
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



# artefacts file comes from Gulay


if h == '0h':
    print ('not implemented')
    sys.exit(1)
elif '13' in h:
    art_cat_file = path+'artefacts_13h_checked.fits'
    art_cat = Table.read(art_cat_file, format='fits')   
    art_cat.rename_column('col0','Source_Name')
    art_cat.add_column(Column(name='artefact_flag',data=np.ones(len(art_cat),dtype=bool)))

    
lofarcat = Table.read(lofarcat_file_srt)

if 'artefact_flag' not in lofarcat.colnames:
    #lofarcat.remove_column('artefact_flag')
    lofarcat.add_column(Column(name='artefact_flag', data=np.zeros(len(lofarcat),dtype=bool)))

lofarcat.sort('Source_Name')
##lofarcat.keep_columns(['Source_Name'])
#tt=join(lofarcat, art_cat, join_type='left', keys=['Source_Name'])
#tt['artefact_flag'].fill_value = 0
#tt = tt.filled()
#tt.sort('Source_Name')


for i in range(len(lofarcat)):
    lofarcat['Source_Name'][i] =  lofarcat['Source_Name'][i].strip()

for s in art_cat['Source_Name']:
    lofarcat['artefact_flag'][lofarcat['Source_Name'] == s] = True

print(np.sum(lofarcat['artefact_flag']),'artefacts')

## write output file
lofarcat.write(lofarcat_file_srt, overwrite=True, serialize_meta=True)
