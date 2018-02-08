#!/usr/bin/python

from astropy.table import Table
import numpy as np

lines=open('double_ID_list.txt').readlines()
ct=Table.read('../LOFAR_HBA_T1_DR1_catalog_v0.99.srl.gmasked.artefacts.fits')
ct.add_index('Source_Name')
st=Table.read('../lgz_v2/LOFAR_HBA_T1_DR1_merge_ID_optical_v0.8.fits')
st.add_index('Source_Name')
ot=st[0:0]
i=0

while i<len(lines):
    l=lines[i].rstrip()
    assert('Group' in l)
    g=l
    i+=1
    sources={}
    while i<len(lines) and lines[i].rstrip()!="":
        bits=lines[i].rstrip().split()
        source=bits[0]
        flag=int(bits[1])
        sources[source]=flag
        i+=1
    i+=1

    # now we have a source list for a group
    s=set([sources[s] for s in sources])
    if s==set((1,311)) or s==set((311,)) or s==set((312,1)) or s==set((600,1)) or s==set((311,600)):
        # pick a source that's in the component catalogue
        j=0
        for s in sources:
            try:
                rc=ct.loc[s]
                r=st.loc[s]
                ot.add_row(r)
                j+=1
            except:
                pass
        if j==0:
            print g,'tzi failed'
            print sources
    else:
        print g,'Not for tzi'
        print sources
    
ot.write('doubles-zoom.fits',overwrite=True)
