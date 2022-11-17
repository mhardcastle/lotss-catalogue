from __future__ import print_function
import numpy as np
import os
import sys
from astropy.table import Table,vstack

# simplified from check_overlap.py since no duplicated sky area

t1=Table.read('Spring/final-v1.0_photoz_v0.2_joined-physical.fits')
t2=Table.read('Fall/final-v1.1_photoz_v0.2_south-physical.fits')

combined=vstack([t1,t2])
combined.sort('RA')
combined.write('combined-release-v0.5.fits',overwrite=True)
