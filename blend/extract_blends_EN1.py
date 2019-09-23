from astropy.table import Table
import numpy as np

t=Table.read('/beegfs/lofar/deepfields/ELAIS_N1_LR/new_optcat_matches/EN1_ML_RUN_fin_overlap_srl_workflow_th.fits')

t=t[t['FLAG_WORKFLOW']==5]

t.write('blends.fits',overwrite=True)
