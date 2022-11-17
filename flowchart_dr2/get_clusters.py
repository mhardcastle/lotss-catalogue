
import sys
import os
import numpy as np
from astropy.table import Table, Column, vstack



if len(sys.argv) == 1:
    print("Usage is : python get_intersecting_sources.py field_code ")
    print('E.g.: python get_intersecting_sources.py 0 ')
    sys.exit(1)

h = str(sys.argv[1])
if 'h' not in h:
    h+='h'
if h not in  ['0h','13h','n0h','n13h','s0h','s13h']:
    print('unknown field code (should be 0h or 13h)',h)
    sys.exit(1)

path = '/Users/w.williams/projects/lofar_surveys/DR2/'
lofarcat_file_srt = path+'LoTSS_DR2_{version}.srl_{h}.lr-full.presort.fits'.format(version=version,h=h)

cat = Table.read(lofarcat_file_srt)

R = 45
out='test_match_c{r:d}.fits'.format(r=R)

if not os.path.isfile(out):
    cmd = "stilts tskymatch2 in1={in1} in2={in2} out={out} join=all1 find=all error=45".format(in1=lofarcat_file_srt,in2=lofarcat_file_srt,out=out)
    os.system(cmd)
    
mcat = Table.read(out)
mcat.keep_columns('Source_Name','GroupSize')

u, ind = np.unique(mcat['Source_Name'], return_index=True)
mcat = mcat[ind]

