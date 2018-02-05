# Process sources with 'host-broken' flag

from astropy.table import Table
from separation import separation
import numpy as np

t=Table.read('HETDEX-LGZ-cat-v0.6-filtered.fits')
t2m=Table.read('/data/lofar/mjh/hetdex_v4/2MASX_hetdex_fix.fits') # Wendy's table
sep2m=60 # arcsec

count=0
problems=np.zeros(len(t),dtype=bool)
for i in range(len(t)):
    r=t[i]
    if r['Hostbroken_prob']<0.4:
        continue
    print r['Source_Name']
    count+=1
    if np.isnan(r['optRA']):
        # host broken up but no host. Hmm...
        print 'No host!'
        sep=3600.0*separation(r['RA'],r['DEC'],t2m['ra'],t2m['dec'])
        j=np.argmin(sep)
        if sep[j]>sep2m:
            print 'No 2MASS source within',sep2m,'arcsec'
            problems[i]=True
            continue
    else:
        sep=3600.0*separation(r['optRA'],r['optDec'],t2m['ra'],t2m['dec'])
        j=np.argmin(sep)
        if sep[j]>sep2m:
            print 'No 2MASS source within',sep2m,'arcsec'
            problems[i]=True
            continue

    print '2MASS source found:',t2m[j]['designation']
    r['optRA']=t2m[j]['ra']
    r['optDec']=t2m[j]['dec']
    r['OptID_Name']='2MASS '+t2m[j]['designation'].rstrip()
    r['Hostbroken_prob']=0
    t[i]=r

print count,'total hostbroken'
print np.sum(problems),'problems'

tp=t[problems]
tp.write('../zoom_v2/hostbroken_problems.fits',overwrite=True)
t.write('HETDEX-LGZ-cat-v0.6-filtered-unbroken.fits',overwrite=True)
