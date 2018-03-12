#!/usr/bin/python

# This is a wrapper round the zoom code which will zoom on a single
# source in a named table. It makes the sub-table with the source in,
# finds the maps, and then calls zoom.py itself.

# paths must be set up as in lgz_create/readme

import os
import sys
from astropy.table import Table
import numpy as np

sourcename=sys.argv[1]
table=sys.argv[2]

t=Table.read(table)
#iterate over table in case of trailing spaces

for i in range(len(t)):
    if t[i]['Source_Name'].rstrip()==sourcename:
        break
else:
    raise RuntimeError('Source not found')

mask=np.zeros(len(t),dtype=bool)
mask[i]=True
t=t[mask]
t.write(sourcename+'-table.fits',overwrite=True)

os.system('python '+os.environ['LGZPATH']+'/utils/make_image_list.py '+sourcename+'-table.fits')
os.system('python '+os.environ['LGZPATH']+'/zoom/zoom.py '+sourcename+'-table.fits flag')
