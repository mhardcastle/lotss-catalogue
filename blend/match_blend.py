# Make a blend input file from noid sources that need to go back to the blend workflow.

from astropy.table import Table
import numpy as np

t=Table.read('sources-v0.2.fits')
ct=Table.read('/beegfs/lofar/deepfields/ELAIS_N1_LR/new_optcat_matches/EN1_ML_RUN_fin_overlap_srl_workflow_th.fits')

filt=t['NoID']==7
filt|=t['NoID']==9

t=t[filt]
print len(t)
filt=np.array([False]*len(ct))
for r in t:
    filt|=ct['Source_Name']==r['Source_Name']

ct[filt].write('blend_noid.fits',overwrite=True)

