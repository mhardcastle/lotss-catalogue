import sys
import numpy as np
from astropy.table import Table


catname=sys.argv[1]
cat = Table.read(catname)
fcflg = sys.argv[2]

print(len(cat),'sources in catalogue',catname)
print(fcflg,'count')
t,i = np.unique(cat[fcflg], return_counts=True)
for tt,ii in zip(t,i): print(tt,ii)
