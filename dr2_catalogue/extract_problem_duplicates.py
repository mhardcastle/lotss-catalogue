from __future__ import print_function
import cPickle as pickle
from astropy.table import Table
import numpy as np
from separation import separation
import os
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors

bits=os.getcwd().split('/')
table=bits[-1].replace('-','_')
print('Table is',table)
con=mdb.connect('127.0.0.1', 'tzi_user', 'IK34daKG', 'duplicates', cursorclass=mdbcursors.DictCursor, autocommit=True)
cur = con.cursor()

zn=pickle.load(open('structure-v0.1-zoomneeded.pickle')) 
zr=pickle.load(open('structure-v0.1-zoomreasons.pickle'))

dups=[]
for n,r in zip(zn,zr):                                   
    if r=='Duplicate ID non-zoom source':
        dups.append(n)

print(len(dups))
        
t=Table.read('sources-v0.1.fits')
filt=np.zeros_like(t,dtype=bool)

for d in dups:
    filt|=t['Source_Name']==d

dt=t[filt]
print(len(dt))

dt.write('duplicates.fits',overwrite=True)
keynames=set([])
for r in dt:
    dist=separation(r['optRA'],r['optDec'],dt['optRA'],dt['optDec'])
    st=dt[dist<1.5/3600.0]
    print(r['Source_Name'],len(st))
    # if there are two (only) sources in st, then we have a pair that we can send to automatic processing. First source will be the lead source and we will store the name of the paired source in the table.
    if len(st)==2 and st[0]['Source_Name'] not in keynames:
        cur.execute('insert into %s values ("%s","%s",NULL)' % (table,st[0]['Source_Name'],st[1]['Source_Name']))
        keynames.add(st[0]['Source_Name'])

con.close()
