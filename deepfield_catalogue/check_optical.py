from astropy.table import Table
import numpy as np
import glob

filename=sorted(glob.glob('final-v*.fits'))[-1]

print 'filename is',filename

t=Table.read(filename)

print 'Total sources',len(t)
print 'Sources with no optRA:',np.sum(np.isnan(t['optRA']))
t[np.isnan(t['optRA'])].write('noid.fits',overwrite=True)
if 'NoID' in t.colnames:
    filter=np.isnan(t['optRA']) & (t['NoID']==0)
    print 'Sources with no optRA and NoID unset:',np.sum(filter)
    t[filter].write('noid2.fits',overwrite=True)
    
print 'Sources with no match:',np.sum(t['NUMBER'].mask)
print 'Sources with optRA but no match:',np.sum(t['NUMBER'].mask)-np.sum(np.isnan(t['optRA']))
print 'Sources with optRA, no match and NoID unset:',np.sum(t['NUMBER'].mask & ~np.isnan(t['optRA']) & (t['NoID']==0))
