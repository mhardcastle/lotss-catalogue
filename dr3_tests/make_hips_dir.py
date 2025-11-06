import glob
import os
import sys

filename=sys.argv[1]

g=glob.glob('/beegfs/car/mjh/healpix-mosaics/*/'+filename)

for f in g[10]:
    bits=f.split('/')
    pixname=bits[-2]
    print('ln -s '+f+' '+pixname+'.fits')
    
