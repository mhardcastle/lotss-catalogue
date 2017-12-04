#!/usr/bin/python

'''
get_visual_flags
add flags to catalogue based on visual inspection of subclasses of sources
'''

from lofar_source_sorter import Mask, Masks_disjoint_complete
import numpy as np
from astropy.table import Table, Column, join
import astropy.coordinates as ac
import astropy.units as u
import os

#################################################################################

path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'
lofarcat_file_srt = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.srl.fixed.presort.fits'



lofarcat = Table.read(lofarcat_file_srt)


#################################################################################
# artefacts come from find_artefact_candidates.py

'''
Candidate artefacts selected by Gulay: TBC - what was the exact selection?
'''
artefactlistfile = 'gg_artefact_case1_3-fixed.fits'
visually_confirmed = False

'''
After doing visual inspection of these candidates: provides new catalogues
'''
artefactlistfile = 'gg_artefact_case1_3-fixed-confirmed.fits'
visually_confirmed = True



artefactlist = Table.read(artefactlistfile)

#select only confirmed ones
if visually_confirmed:
    artefactlist = artefactlist[artefactlist['visual_flag']==1]

# for now, no artefacts
artefact = np.zeros(len(lofarcat),dtype=bool)
if 'artefact_flag' not in lofarcat.colnames:
    lofarcat.add_column(Column(artefact,'artefact_flag'))
else:
    #rewrite artefact info
    lofarcat['artefact_flag'] *= False
for n in artefactlist['Source_Name']:
    ni = np.where(lofarcat['Source_Name']==n)[0][0]
    lofarcat['artefact_flag'][ni] = True    


#################################################################################
# the following come from outputs of classify.py on various subsamples

#################################################################################
### nhuge_2masx
#1: bright galaxy
#2: complex
#3: no match
#4: artefact

nhuge_2masx_vc_cat_file = 'fullsample/sample_all_src_clean_large_faint_nhuge_2masx-vflag.fits'
nhuge_2masx_vc_cat = Table.read(nhuge_2masx_vc_cat_file)

if 'nhuge_2masx_flag' in lofarcat.colnames:
    lofarcat.remove_column('nhuge_2masx_flag')
lofarcat.sort('Source_Name')
tt=join(lofarcat, nhuge_2masx_vc_cat, join_type='left')
tt['visual_flag'].fill_value = 0
tt = tt.filled()
tt.sort('Source_Name')
tt.rename_column('visual_flag','nhuge_2masx_flag')


lofarcat.add_column(tt['nhuge_2masx_flag'])

#################################################################################
### clustered
#1: artefact
#2: complex

clustered_vc_cat_file = 'fullsample/sample_all_src_clean_small_nisol_clustered-vflag.fits'
clustered_vc_cat = Table.read(clustered_vc_cat_file)

if 'clustered_flag' in lofarcat.colnames:
    lofarcat.remove_column('clustered_flag')
lofarcat.sort('Source_Name')
tt=join(lofarcat, clustered_vc_cat, join_type='left')
tt['visual_flag'].fill_value = 0
tt = tt.filled()
tt.sort('Source_Name')
tt.rename_column('visual_flag','clustered_flag')


lofarcat.add_column(tt['clustered_flag'])


#################################################################################
### large faint clustered
#1: artefact
#2: complex

Lclustered_vc_cat_file = 'testsample_large/sample_all_src_clean_large_faint_nhuge_n2masx_nisol_clustered-vflag.fits'
Lclustered_vc_cat = Table.read(Lclustered_vc_cat_file)

if 'Lclustered_flag' in lofarcat.colnames:
    lofarcat.remove_column('Lclustered_flag')
lofarcat.sort('Source_Name')
tt=join(lofarcat, Lclustered_vc_cat, join_type='left')
tt['visual_flag'].fill_value = 0
tt = tt.filled()
tt.sort('Source_Name')
tt.rename_column('visual_flag','Lclustered_flag')


tt['Lclustered_flag'][tt['Source_Name']=='ILTJ135637.88+473205.2'] = 2

lofarcat.add_column(tt['Lclustered_flag'])


#################################################################################
### huge faint
#1: send to LGZ
#2: bright galaxy match
#3: no prospect of ID
#4: artefact

huge_faint_vc_cat_file = 'fullsample/sample_all_src_clean_large_faint_huge-vflag.fits'
huge_faint_vc_cat = Table.read(huge_faint_vc_cat_file)

if 'huge_faint_flag' in lofarcat.colnames:
    lofarcat.remove_column('huge_faint_flag')
lofarcat.sort('Source_Name')
tt=join(lofarcat, huge_faint_vc_cat, join_type='left')
tt['visual_flag'].fill_value = 0
tt = tt.filled()
tt.sort('Source_Name')
tt.rename_column('visual_flag','huge_faint_flag')


## fix a few bad ones!!
tt['huge_faint_flag'][tt['Source_Name']=='ILTJ141956.60+533054.4'] = 3
tt['huge_faint_flag'][tt['Source_Name']=='ILTJ123133.59+484958.6'] = 3
tt['huge_faint_flag'][tt['Source_Name']=='ILTJ135431.79+542009.6'] = 1
tt['huge_faint_flag'][tt['Source_Name']=='ILTJ114328.67+524240.1'] = 1
tt['huge_faint_flag'][tt['Source_Name']=='ILTJ130840.28+540437.0'] = 1
tt['huge_faint_flag'][tt['Source_Name']=='ILTJ132513.22+535113.5'] = 1
tt['huge_faint_flag'][tt['Source_Name']=='ILTJ133233.91+541927.4'] = 1
tt['huge_faint_flag'][tt['Source_Name']=='ILTJ132919.17+530505.1'] = 1
tt['huge_faint_flag'][tt['Source_Name']=='ILTJ105949.84+534811.6'] = 1

lofarcat.add_column(tt['huge_faint_flag'])


#################################################################################
### large, not huge faint
#(1) Send to LGZ
#(2) Accept ML match
#(3) No good match
#(4) Too zoomed in
#(5) Artefact

huge_faint_vc_cat_file = 'toclassify_171124//large_faint_toclassify-vflag.fits'
huge_faint_vc_cat = Table.read(huge_faint_vc_cat_file)

if 'nhuge_faint_flag' in lofarcat.colnames:
    lofarcat.remove_column('nhuge_faint_flag')
lofarcat.sort('Source_Name')
tt=join(lofarcat, huge_faint_vc_cat, join_type='left')
tt['visual_flag'].fill_value = 0
tt = tt.filled()
tt.sort('Source_Name')
tt.rename_column('visual_flag','nhuge_faint_flag')


lofarcat.add_column(tt['nhuge_faint_flag'])


#################################################################################


## write output file

if os.path.exists(lofarcat_file_srt):
    os.remove(lofarcat_file_srt)
lofarcat.write(lofarcat_file_srt)