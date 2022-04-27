import glob
import os
from astropy.table import Table
import sys
import numpy as np

optcat='/beegfs/lofar/mjh/dr2/dr2_combined.fits'
infile='mksp-output.fits' # with original names
mergeout='mksp-merged.fits'

os.system('/soft/topcat/stilts tskymatch2 ra1=optRA dec1=optDec ra2=RA dec2=DEC error=3.0 join=all1 in1=%s in2=%s out=%s' % (infile,optcat,mergeout))

t=Table.read(mergeout)
t['RA_1'].name='RA'
t['DEC_1'].name='DEC'
t['RA_2'].name='idRA'
t['DEC_2'].name='idDEC'

t.write('mksp-fixed.fits',overwrite=True)
