from __future__ import print_function
from astropy.table import Table
import glob
import os
import numpy as np
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors

# Merge prefilter for the Spring field from the various databases --
# this avoids needing to ask Wendy to do it manually

ct=Table.read('source_lr.fits')
sd={r['Source_Name']:i for i,r in enumerate(ct)}

con=mdb.connect('127.0.0.1', 'prefilter_user', 'WQ98xePI', 'prefilter', cursorclass=mdbcursors.DictCursor, autocommit=True)
cur = con.cursor()

for table in ['13h40','13h60','13h60b','13h60fix','8h','hetdex','last','mfix','missing','missing2']:
    print('=== %s ===' %table)
    cur.execute('select * from '+table)
    results=list(cur.fetchall())
    for r in results:
        source=r['object']
        if source in sd:
            print('Setting',source,'to',r['classification'])
            ct[sd[source]]['Prefilter']=r['classification']
        else:
            print(source,'not found!')

ct['Blend_prefilter']=0

for table in ['blend_filter_Spring','blend_filter_Fall']:
    print('=== %s ===' %table)
    cur.execute('select * from '+table)
    results=list(cur.fetchall())
    for r in results:
        source=r['object']
        if source in sd:
            print('Setting',source,'to',r['classification'])
            ct[sd[source]]['Blend_prefilter']=r['classification']
        else:
            pass
            #print(source,'not found!')
    
            
ct.write('source_lr_prefilter.fits',overwrite=True)
