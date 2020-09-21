from __future__ import print_function
from astropy.table import Table
import glob
import os

from scale_fits import scale_fits

wd=os.getcwd()
if 'Bootes' in wd:
    scale=1.0/0.859
elif 'Lockman' in wd:
    scale=1.0/0.920
else:
    raise NotImplementedError('Cannot scale this field, what is it?')

g=glob.glob('*.fits')

for f in g:
    if '.cat.' in f or 'pybdsmmask' in f or 'beam' in f:
        print('Skipping',f)
    else:
        if 'scaled' not in f:
            outname=f.replace('.fits','.scaled.fits')
        else:
            outname=f.replace('scaled','rescaled')
        print('Scaling',f,'to',outname)

        scale_fits(f,scale,outname)

for f in g:
    if '.cat.' in f:
        t=Table.read(f)
        outname=f.replace('scaled','rescaled')
        for c in t.colnames:
            if 'flux' in c or 'rms' in c or 'mean' in c:
                t[c]/=scale
        t.write(outname)
        
        
