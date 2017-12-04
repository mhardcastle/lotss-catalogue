#!/usr/bin/python

'''
find_artefacts
flag artefacts based on candidate selection, visual inspection (requires catalogues provided through other means)
'''

from lofar_source_sorter import Mask, Masks_disjoint_complete
import numpy as np
from matplotlib import patches
from astropy.table import Table, Column
import astropy.coordinates as ac
import astropy.units as u
import os

#################################################################################

path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'
lofarcat_file_srt = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.srl.fixed.presort.fits'


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


#################################################################################

lofarcat = Table.read(lofarcat_file_srt)

artefactlist = Table.read(artefactlistfile)

#select only confirmed ones
if visually_confirmed:
    artefactlist = artefactlist[artefactlist['visual_flag']==1]

# for now, no artefacts
artefact = np.zeros(len(lofarcat),dtype=bool)
if 'artefact' not in lofarcat.colnames:
    lofarcat.add_column(Column(artefact,'artefact'))
else:
    #rewrite artefact info
    lofarcat['artefact'] *= False
for n in artefactlist['Source_Name']:
    ni = np.where(lofarcat['Source_Name']==n)[0][0]
    lofarcat['artefact'][ni] = True    


## write output file

if os.path.exists(lofarcat_file_srt):
    os.remove(lofarcat_file_srt)
lofarcat.write(lofarcat_file_srt)