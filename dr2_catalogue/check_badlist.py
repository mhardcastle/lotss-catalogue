from __future__ import print_function
import pickle
import glob
from astropy.table import Table

bad=pickle.load(open('badlist.pickle'))

bad=sorted(list(set(bad)))
g=sorted(glob.glob('sources-v*.fits'))

filename=g[-1]
print('loading',filename)
t=Table.read(filename)

for s in bad:
    r=t[t['Source_Name']==s][0]
    print(s,r['Created'],r['Position_from'],r['optRA'])
    
