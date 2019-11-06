from astropy.table import Table
import numpy as np
from separation import separation
import os
import glob

dir=os.getcwd()
field=os.path.basename(dir)
print 'field is',field

if field=='bootes':
    ct=Table.read('/beegfs/lofar/deepfields/Bootes_LR/new_fdeep_matches/Bootes_ML_RUN_fin_overlap_srl_workflow_th.fits')
elif field=='lockman':
    ct=Table.read('/beegfs/lofar/deepfields/Lockman_LR/updated_LR_cols/LH_ML_RUN_fin_overlap_srl_workflow_th.fits')
elif field=='en1':
    ct=Table.read('/beegfs/lofar/deepfields/ELAIS_N1_LR/new_optcat_matches/EN1_ML_RUN_fin_overlap_srl_workflow_th.fits')
else:
    print 'Not in correct working directory'
    sys.exit(1)


g=sorted(glob.glob('final-v*.fits'))
infile=g[-1]
print 'Processing',infile

t=Table.read(infile)
seps=[]

for r in t:
    sep=3600.0*separation(r['RA'],r['DEC'],r['ALPHA_J2000'],r['DELTA_J2000'])
    seps.append(sep)

t['sep']=seps
t.write(infile.replace('.fits','-withsep.fits'),overwrite=True)

filter=t['FLAG_WORKFLOW']==1
filter&=t['sep']>3

nt=t[filter]

filter=np.array([False]*len(ct))
for r in nt:
    filter|=(ct['Source_Name']==r['Source_Name'])

ct[filter].write('bigsep.fits',overwrite=True)

filter=t['NoID']==7
filter|=t['NoID']==9

nt=t[filter]
filter=np.array([False]*len(ct))
for r in nt:
    filter|=(ct['Source_Name']==r['Source_Name'])

ct[filter].write('noid.fits',overwrite=True)
