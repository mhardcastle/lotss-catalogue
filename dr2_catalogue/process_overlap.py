#!/usr/bin/env python

from __future__ import print_function
import os
from astropy.table import Table
import cPickle as pickle
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from separation import separation
from tqdm import tqdm
from astropy_healpix import HEALPix
from astropy import units as u
from multiprocessing import Pool, cpu_count

def select_multi(i):
    if np.isnan(optras[i]): return []
    thishp=hd[hps[i]] # the ones with matching healpix
    if len(thishp)==1: return []
    dist=separation(optras[i],optdecs[i],thishp[1],thishp[2])
    filt=dist<1.5/3600.0
    sm=np.sum(filt)
    if sm>1:
        return list(thishp[0][filt])
    else:
        return []

def slice_multi(hpv):
    filt=(hps==hpv)
    return (names[filt],optras[filt],optdecs[filt])

print('*** process_overlap starting ***')
with open('optical.pickle') as pf:
    (names,optras,optdecs)=pickle.load(pf)

hd={}
hp = HEALPix(nside=512)
hps=hp.lonlat_to_healpix(optras*u.deg,optdecs*u.deg)

p=Pool(cpu_count())

hpl=list(set(hps))
for i,result in enumerate(tqdm(p.imap(slice_multi,hpl),total=len(hpl),desc='Make dictionary')):
    hd[hpl[i]]=result
p.close()
del(p)
    
p=Pool(48)
badlist=[]
for result in tqdm(p.imap(select_multi,range(len(optras))),total=len(optras),desc='Find neighbours'):
    badlist+=result

p.close()
del(p)
print(len(badlist),'duplicates found')

with open('badlist.pickle','w') as pf:
    pickle.dump(badlist,pf)
    
