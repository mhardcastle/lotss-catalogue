# Find a source name in the catalogues and figure out how to do the
# appropriate zoom

import sys
from astropy.table import Table
import glob

name=sys.argv[1]

t=Table.read('/beegfs/general/lofar/LOFAR_HBA_T1_DR1_merge_ID_optical_v1.0.fits')

tc=Table.read('/beegfs/general/lofar/LOFAR_HBA_T1_DR1_merge_ID_v1.0.comp.fits')

s=t[t['Source_Name']==name]
if len(s)!=1:
    raise RuntimeError('Source is not in main catalogue')

r=s[0]
print 'Source',name,'has ID_flag',r['ID_flag']

ct=tc[tc['Source_Name']==name]

if len(ct)==0:
    raise RuntimeError('Source has no components!')
elif len(ct)==1:
    print name,'is a single-component source'
else:
    print name,'has',len(ct),'components'
    for c in ct:
        print c['Component_Name']

# Now try to find out whether it is in LGZ/zoom already

for d in ['/data/lofar/mjh/hetdex_v4/zoom','/data/lofar/mjh/hetdex_v4/zoom_v2']:
    lines=open(d+'/components.txt').readlines()
    for l in lines:
        bits=l.rstrip().split()
        if name==bits[0]:
            print 'Source name appears as composite source name in',d
        for c in ct:
            if c['Component_Name']==bits[1]:
                print 'Component name appears as component of',bits[0],'in',d

# Now check the zoom files

for d in ['/data/lofar/mjh/hetdex_v4/zoom','/data/lofar/mjh/hetdex_v4/zoom_v2']:
    g=glob.glob(d+'/*.txt')
    for f in g:
        if 'list' in f or 'components.txt' in f:
            continue
        lines=open(f).readlines()
        for l in lines:
            for c in ct:
                if c['Component_Name'] in l:
                    print 'Component appears in zoom file',f

