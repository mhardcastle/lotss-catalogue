from __future__ import print_function
import sys
from astropy.table import Table
import os
import numpy as np

infile=open('radio-galaxy-zoo-lofar-classifications.csv')


ramin=326
ramax=60
decmin=16
decmax=36

delta=0.1
deltara=delta*np.cos(np.pi*(decmin+decmax)/360.0)


outname='Fall'
#lrcat='/beegfs/lofar/wwilliams/lofar_surveys/DR2/lr/LoTSS_DR2_v100.srl_13h.lr-full.fits'
lrcat='/beegfs/lofar/wwilliams/lofar_surveys/DR2/flowchart_dec7/LoTSS_DR2_v100.srl_0h.lr-full.sorted_step3_flux4.fits'
gaucat='/beegfs/lofar/wwilliams/lofar_surveys/DR2/lr/LoTSS_DR2_v100.gaus_0h.lr-full.fits'
if not os.path.isdir(outname):
    os.system('mkdir '+outname)
os.chdir(outname)
print('Filter RGZ')
outfile=open(outname+'.csv','w')

for i,line in enumerate(infile):

    if (i % 10000)==0:
        print('.',end='')
        sys.stdout.flush()
    if i==0:
        outfile.write(line)
    else:
        bits=line.split(',')
        if bits[4]!="11973":
            continue
        # now locate RA and DEC
        ra=None
        dec=None
        for b in bits[5:]:  
            if '""ra""' in b:
                ra=float(b[7:].replace('"',''))
            if '""dec""' in b:
                dec=float(b[8:].replace('"',''))
        if ra is None or dec is None:
            continue
            #print bits
            #raise RuntimeError("Failed to find the position")
        #print bits
        #print ra,ramin,ramax,dec,decmin,decmax
        #stop
        if (ra>=ramin or ra<=ramax) and dec>=decmin and dec<=decmax:
            outfile.write(line)

# apply the same cuts to the optical catalogue
# this will eventually be one (or more) of the augmented optical catalogues

outfile.close()
print()
print('Filter optical')
optcat='/beegfs/lofar/mjh/dr2/dr2_combined.fits' # eventually annotated version
t=Table.read(optcat)
filt=t['RA']>=ramin-deltara
filt|=t['RA']<=ramax+deltara
filt&=t['DEC']>=decmin-delta
filt&=t['DEC']<=decmax+delta
t[filt].write('optical.fits',overwrite=True)

# apply the cuts to the LR catalogues
print('Filter LR')
t=Table.read(lrcat)
filt=t['RA']>=ramin-deltara
filt|=t['RA']<=ramax+deltara
filt&=t['DEC']>=decmin-delta
filt&=t['DEC']<=decmax+delta
t[filt].write('source_lr.fits',overwrite=True)

print('Filter Gaussian')
t=Table.read(gaucat)
filt=t['RA']>=ramin-deltara
filt|=t['RA']<=ramax+deltara
filt&=t['DEC']>=decmin-delta
filt&=t['DEC']<=decmax+delta
t[filt].write('gaussian_lr.fits',overwrite=True)

# We allow a border because associations may pick up sources outside the cut area. We'll then have to cut down the output catalogue by optical position at the end of the process.
