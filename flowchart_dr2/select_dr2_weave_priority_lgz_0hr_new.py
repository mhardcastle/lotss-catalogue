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

from lofar_source_sorter_dr2 import Mask, Masks_disjoint_complete


#################################################################################


priority = '1'
selno = '2'


path = '/data2/wwilliams/projects/lofar_surveys/LoTSS-DR2-Feb2020/'
#lofarcat_file_srt = path+'LoTSS_DR2_v100.srl_{h}.lr-full.sorted_step1.fits'.format(h=h)
lofarcat_file_srt = path+'LoTSS_DR2_v100.srl_0h.lr-full.sorted_step2_flux4.hdf5'
sel_file = path+'lgz_selection_nov18/LoTSS_DR2_v100.srl_0h.lr-full.sorted_step2_flux4.lgz_selection_{i}.fits'.format(i=selno)


lofarcat = Table.read(lofarcat_file_srt)


sel_pri = (lofarcat['WEAVE_priority{i}'.format(i=priority)] == True)
#sel = (lofarcat['WEAVE_priority1'] == True) & (lofarcat['WEAVE_priority1a'] == False)  # special case of 1 and not 1a  make up 1b


# flowchart now accounts for this
#sel =  sel_pri & (
        ##(lofarcat['FC_flag2'] == 5) |    # large  4<S<8 mJy, ML=LGZ  ## NB already done
        #((lofarcat['FC_flag2'] == 12) & (lofarcat['ML_flag']==1)) |   # isol nS S>8mJy msource=LGZ  -## NB all ML=LGZ already done, so add the LR ones
        ##(lofarcat['FC_flag2'] == 15) |   # n-isol, clustered, S>4mJy, ML=LGZ   ## NB clustered already done with LGZ
        #(lofarcat['FC_flag2'] == 26)  |   # n-isol, nclustered, S, frat, sep, S>4mJy, ML=LGZ
        #((lofarcat['FC_flag2'] == 20)& (lofarcat['ML_flag']==1)) )    # n-isol, nclustered, nS, S>4mJy, msource=LGZ  ## NB all ML=LGZ already done, so add the LR ones

sel =  sel_pri & (
        (lofarcat['FC_flag2'] == 6) | \
        (lofarcat['FC_flag2'] == 14) | \
        (lofarcat['FC_flag2'] == 21) | \
        (lofarcat['FC_flag2'] == 27))


lofarcat1 = lofarcat[sel]

lofarcat1.write(sel_file, overwrite=True)

sel_file = path+'lgz_selection_nov18/LoTSS_DR2_v100.srl_0h.lr-full.sorted_step2_flux4.prefilter_lgz_selection_{i}.fits'.format(i=selno)

#sel = sel_pri & (
        #(lofarcat['FC_flag2'] == 6) |    # large  4<S<8 mJy, ML=LR
        #((lofarcat['FC_flag2'] == 13) & (lofarcat['ML_flag']==1)) |   # isolate S S>8mJy msource=prefilt   ## NB ML=LGZ already done in lgz
        ##(lofarcat['FC_flag2'] == 16) |   # n-isol, clustered, S>4mJy, ML=LR   ## NB clustered already done with LGZ
        #(lofarcat['FC_flag2'] == 27) |    # n-isol, nclustered, S, frat, sep, S>4mJy, ML=LR
        #((lofarcat['FC_flag2'] == 21)  & (lofarcat['ML_flag']==1)) )   # n-isol, nclustered, nS, S>4mJy, msource=prefilt   ## NB ML=LGZ already done in lgz



sel =  sel_pri & (
        (lofarcat['FC_flag2'] == 5) | \
        (lofarcat['FC_flag2'] == 9) | \
        (lofarcat['FC_flag2'] == 13) | \
        (lofarcat['FC_flag2'] == 15) | \
        (lofarcat['FC_flag2'] == 16) | \
        (lofarcat['FC_flag2'] == 16) )



lofarcat1 = lofarcat[sel]

lofarcat1.write(sel_file, overwrite=True)

'''
Box 7: compact, clustered
1536 sources
select with (t['FC_flag3']==7) 
-- all clustered went to LGZ  * q? what were the outcomes?
-- no prefilter, but would only be 20 sources

Box 3: S< 8mJy and size > 15" (47,834 srcs) - some details on the numbers here:
 14,154 : 33,680  (S>4 : <= 4 mJy)
13,125: 1,026   :  22,576:11,104  (ML LGZ : ML LR) (ML=machine learning)
so this would select 13,125 sources now for LGZ (4<S<8 mJy, size>15" and ML=LGZ)
select with (t['FC_flag3']==3) & (t['Total_flux']>4) & (t['ML_flag']==0)
- Note that (4,039 of these are >30" so have no LR information)
* 1029 for prefilter
FC2_flag2==6 

Msources: with ML=LGZ
(t['FC_flag3']==6) &  (t['ML_flag']==0)
(t['FC_flag3']==8) &  (t['ML_flag']==0)
6,244 and 9,972 srcs

* now isol and >8mJy 
-> prefilter 10363 (8592:1771  LR:LGZ)    +++ 8592 for prefilter
FC2_flag2==13 & ML_flag==1
-> lgz 404 (274:130)   +++ 274 for LGZ directly
FC2_flag2==12 & ML_flag==1

* now non-isol and >4 mJy 
-> prefilter 6706 (2387:4319 LR:LGZ)   +++ 2387 for prefilter
FC2_flag2==21 & ML_flag==1
-> lgz 286 (46:240)   +++ 46 for LGZ directly
FC2_flag2==20 & ML_flag==1



* tail end 
-> 1339 for prefilter
FC2_flag2==27
-> 3176 for LGZ
FC2_flag2==26

** also special case - do not need to do any msource and clustered and > 4mJy by prefilter/lgz because they were already done by lgz   ### actually only 2 already in lgz ... so ignore
'''
