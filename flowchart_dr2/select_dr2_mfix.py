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

#from lofar_source_sorter_dr2 import Mask, Masks_disjoint_complete


#################################################################################


hetdex_fields = ['P10Hetdex', 'P11Hetdex12', 'P12Hetdex11', 'P14Hetdex04', 'P15Hetdex13', 'P164+55', 'P169+55', 'P16Hetdex13', 'P173+55', 'P178+55', 'P182+55', 'P187+55', 'P18Hetdex03', 'P191+55', 'P196+55', 'P19Hetdex17', 'P1Hetdex15', 'P200+55', 'P205+55', 'P206+50', 'P206+52', 'P209+55', 'P21', 'P210+47', 'P211+50', 'P213+47', 'P214+55', 'P217+47', 'P218+55', 'P219+50', 'P219+52', 'P221+47', 'P223+50', 'P223+52', 'P223+55', 'P225+47', 'P227+50', 'P227+53', 'P22Hetdex04', 'P23Hetdex20', 'P25Hetdex09', 'P26Hetdex03', 'P27Hetdex09', 'P29Hetdex19', 'P30Hetdex06', 'P33Hetdex08', 'P34Hetdex06', 'P35Hetdex10', 'P37Hetdex15', 'P38Hetdex07', 'P39Hetdex19', 'P3Hetdex16', 'P41Hetdex', 'P42Hetdex07', 'P4Hetdex16', 'P6', 'P7Hetdex11', 'P8Hetdex']

#priority = '1'
#priority = '2'


path = '/Users/w.williams/projects/lofar_surveys/DR2/'
#lofarcat_file_srt = path+'LoTSS_DR2_{version}.srl_{h}.lr-full.sorted_step1.fits'.format(version=version,h=h)
lofarcat_file_srt = path+'LoTSS_DR2_{version}.srl_13h.lr-full.sorted_step2_flux4.hdf5'


lofarcat = Table.read(lofarcat_file_srt)

lofarcat.add_column(Column(name='hetdex_field',data=np.zeros(len(lofarcat))))
for field in hetdex_fields:
    lofarcat['hetdex_field'] = lofarcat['hetdex_field'] + (lofarcat['Mosaic_ID'] == field)

for priority in ['hetdex','hetdex_missing','last','1','2','3']:
    sel_file = path+'lgz_selection_sep9/LoTSS_DR2_{version}.srl_13h.lr-full.sorted_step2_flux4.lgz_weave_selection_{i}.fits'.format(i=priority)
    if priority in ['1' , '2' , '1a', '3']:
        sel_pri = (lofarcat['WEAVE_priority{i}'.format(i=priority)] == True)
    elif priority == '1b':
        sel_pri = (lofarcat['WEAVE_priority1'] == True) & (lofarcat['WEAVE_priority1a'] == False)  # special case of 1 and not 1a  make up 1b
    # fields not included in dr1 
    elif priority == 'hetdex_missing':
        sel_pri = (lofarcat['Mosaic_ID']=='P210+52') | (lofarcat['Mosaic_ID']=='P214+52') | (lofarcat['Mosaic_ID']=='P215+50') | (lofarcat['Mosaic_ID']=='P31Hetdex19')
    elif priority == 'hetdex':
        sel_pri = (lofarcat['hetdex_field']>0)
    # remaining non-priority areas
    elif priority == 'last':
        sel_pri = (~lofarcat['WEAVE_priority1']) & (~lofarcat['WEAVE_priority2']) & (~lofarcat['WEAVE_priority3']) & (lofarcat['hetdex_field']==0)

    sel =  sel_pri & (
            (lofarcat['FC_flag2'] == 5) | \
            (lofarcat['FC_flag2'] == 12) | \
            (lofarcat['FC_flag2'] == 15) | \
            (lofarcat['FC_flag2'] == 26) | \
            (lofarcat['FC_flag2'] == 20) )



    lofarcat1 = lofarcat[sel]

    lofarcat1.write(sel_file, overwrite=True)

    sel_file = path+'lgz_selection_sep9/LoTSS_DR2_{version}.srl_13h.lr-full.sorted_step2_flux4.prefilter_lgz_weave_selection_{i}.fits'.format(i=priority)

    sel = sel_pri & (
            (lofarcat['FC_flag2'] == 6) | \
            (lofarcat['FC_flag2'] == 13) | \
            (lofarcat['FC_flag2'] == 16) | \
            (lofarcat['FC_flag2'] == 27) | \
            (lofarcat['FC_flag2'] == 21) )



    lofarcat1 = lofarcat[sel]

    lofarcat1.write(sel_file, overwrite=True)
