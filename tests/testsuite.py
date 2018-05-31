# Run all tests

# Call with the names of source and component tables

import sys
from astropy.table import Table, join
import numpy as np
from utils import banner,test_duplicate
import astropy.table
from progressbar import ProgressBar
pbar = ProgressBar()

def test_duplicate_id(t):
    banner('Testing for duplicate optical IDs')
    count=test_duplicate(t,'ID_name','optical ID name',exclude=['Mult','Altered'])
    print 'Found',count,'cases of duplication'
    return count==0

def test_duplicate_sourcename(t):
    banner('Testing for optical cat sourcenames')
    count=test_duplicate(t,'Source_Name','source name')
    print 'Found',count,'cases of duplication'
    return count==0

def test_duplicate_compname(t):
    banner('Testing for duplicate component cat component names')
    count=test_duplicate(t,'Component_Name','component name')
    print 'Found',count,'cases of duplication'
    return count==0

def test_catalogue_mismatch(t,tc):
    banner('Testing for source name mismatches')
    count=0
    print 'Making lists...'
    tn=list(set(t['Source_Name']))
    tcn=list(set(tc['Source_Name']))
    print 'Sorting...'
    tn.sort()
    tcn.sort()
    print 'Checking...'
    i=0
    j=0
    while i<len(tn) and j<len(tcn):
        if tn[i]==tcn[j]:
            i+=1
            j+=1
        elif tn[i]<tcn[j]:
            name=tn[i]
            source=t[t['Source_Name']==name][0]
            print '%s (%i) is in source catalogue but not component catalogue' % (name,source['ID_flag']),i,j
            count+=1
            i+=1
        elif tcn[j]<tn[i]:
            print tcn[j],'is in component catalogue but not source catalogue',i,j
            j+=1
            count+=1
        else:
            raise RuntimeError('Cannot happen')
    print 'Found',count,'mismatches'
    return count==0
    
def test_idflag_mismatch(t,tc):
    banner('Testing for mismatched ID flags')
    tt = join(t,tc, keys=['Source_Name'])
    count = np.sum(tt['ID_flag_1'] != tt['ID_flag_2'])
    print 'Found',count,'cases of mismatching ID_flag\'s'
    return count==0

def test_missing_source(tc,ta,ts):
    banner('Testing for missing sources from the PyBDSF catalogue')
    count=0
    for r in pbar(ts):
        name=r['Source_Name']
        if np.sum(tc['Component_Name']==name)==0:
            # not in component catalogue
            if np.sum(tc['Deblended_from']==name)==0:
                # not deblended
                if np.sum(ta['Source_Name']==name)==0:
                    print 'PyBDSF source',name,'not found'
                    count += 1
    print 'Found',count,'missing sources'
    return count==0
                    
t=Table.read(sys.argv[1])
tc=Table.read(sys.argv[2])
ta=Table.read(sys.argv[3])
ts=Table.read(sys.argv[4])

test_duplicate_id(t)
test_missing_source(tc,ta,ts)
test_duplicate_sourcename(t)
test_duplicate_compname(tc)
test_catalogue_mismatch(t,tc)
test_idflag_mismatch(t,tc)
