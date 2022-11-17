#!/usr/bin/env python3

'''
make the selection of sources for Martin to feed into LGZ
weave priority , mostly flux above 4
'''
import sys
import os
import numpy as np
import glob

from astropy.table import Table, Column, join, vstack
import astropy.coordinates as ac
import astropy.units as u

from lofar_source_sorter_dr2 import Mask, Masks_disjoint_complete


#################################################################################


priority = '1a'
priority = '1b'
#priority = '2'


path = '/data2/wwilliams/projects/lofar_surveys/DR2/'
# run in this dir:
#pathpre = '/data2/wwilliams/projects/lofar_surveys/DR2/lgz_selection'

lofarcat_file_srt = path+'LoTSS_DR2_{version}.srl_13h.lr-full.sorted_step2_flux4.hdf5'

flist = sorted(glob.glob('*13h*selection*fits'))

sadd_lgz = []
sadd_prefilter = []

for f in flist:
    t = Table.read(f)
    tpre = Table.read('prebug/'+f)
    print(f, len(t), len(tpre))
    
    s_not_in_pre = [name for name in t['Source_Name'] if name not in tpre['Source_Name']]
    s_not_in_new = [name for name in tpre['Source_Name'] if name not in t['Source_Name']]
    
    print(s_not_in_pre)
    print(s_not_in_new)
    
    if 'prefilter_lgz' in f:
        sadd_prefilter +=s_not_in_pre
    elif 'lgz' in f:
        sadd_lgz += s_not_in_pre


sel_file = path+'LoTSS_DR2_{version}.srl_0h.lr-full.sorted_step2_flux4.lgz_selection_thresholdfix.fits'
lofarcat = Table.read(lofarcat_file_srt)

for i  in range(len(lofarcat)):
     lofarcat['Source_Name'][i] = lofarcat['Source_Name'][i].strip()

minds = []
for s in sadd_lgz:
    print (s, s in lofarcat['Source_Name'])
    minds.append(np.where(lofarcat['Source_Name']==s)[0][0])
lofarcat1 = lofarcat[minds]

lofarcat1.write(sel_file, overwrite=True)


sel_file = path+'LoTSS_DR2_{version}.srl_0h.lr-full.sorted_step2_flux4.prefilter_selection_thresholdfix.fits'
lofarcat = Table.read(lofarcat_file_srt)
minds = []
for s in sadd_prefilter:
    minds.append(np.where(lofarcat['Source_Name']==s)[0][0])
lofarcat1 = lofarcat[minds]

lofarcat1.write(sel_file, overwrite=True)


#if priority == '1' or priority == '2' or priority == '1a':
    #sel_pri = (lofarcat['WEAVE_priority{i}'.format(i=priority)] == True)
#elif priority == '1b':
    #sel_pri = (lofarcat['WEAVE_priority1'] == True) & (lofarcat['WEAVE_priority1a'] == False)  # special case of 1 and not 1a  make up 1b

#sel =  sel_pri & (
        #(lofarcat['FC_flag2'] == 5) | \
        #(lofarcat['FC_flag2'] == 12) | \
        #(lofarcat['FC_flag2'] == 15) | \
        #(lofarcat['FC_flag2'] == 26) | \
        #(lofarcat['FC_flag2'] == 20) )



#lofarcat1 = lofarcat[sel]

#lofarcat1.write(sel_file, overwrite=True)

#sel_file = path+'LoTSS_DR2_{version}.srl_13h.lr-full.sorted_step2_flux4.prefilter_lgz_weave_selection_{i}.fits'.format(i=priority)

#sel = sel_pri & (
        #(lofarcat['FC_flag2'] == 6) | \
        #(lofarcat['FC_flag2'] == 13) | \
        #(lofarcat['FC_flag2'] == 16) | \
        #(lofarcat['FC_flag2'] == 27) | \
        #(lofarcat['FC_flag2'] == 21) )



#lofarcat1 = lofarcat[sel]

#lofarcat1.write(sel_file, overwrite=True)
