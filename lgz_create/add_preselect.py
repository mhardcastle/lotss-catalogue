# Add the preselect choices to a column in a FITS table and write out again

from astropy.table import Table
import sys
import numpy as np

t=Table.read(sys.argv[1])
outname=sys.argv[2]
t['Prefilter']=0

for f in sys.argv[3:]:
    lines=[l.rstrip() for l in open(f).readlines()]
    for l in lines:
        bits=l.split(',')
        match=t['Source_Name']==bits[1]
        if np.sum(match)!=1:
            raise RuntimeError('Failed to find source '+bits[1])
        i=np.argmax(match)
        print i,bits
        t[i]['Prefilter']=int(bits[0])

t.write(outname,overwrite=True)

           
