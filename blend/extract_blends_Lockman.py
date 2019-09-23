from astropy.table import Table
import numpy as np

t=Table.read('/beegfs/lofar/deepfields/Lockman_LR/updated_LR_cols/LH_ML_RUN_fin_overlap_srl_workflow_th.fits')

t=t[t['FLAG_WORKFLOW']==5]

t.write('blends.fits',overwrite=True)
