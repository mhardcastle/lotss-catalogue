from __future__ import print_function
import glob
from astropy.table import Table
import numpy as np

infile=sorted(glob.glob('sources-v*.fits'))[-1]
print('Reading',infile)
t=Table.read(infile)

filt=np.isnan(t['optRA'])
filt|=t['Position_from']=="Visual inspection"

print('Writing filtered output')
t[filt].write('ridgeline-in.fits',overwrite=True)
