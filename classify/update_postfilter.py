from __future__ import print_function
from astropy.table import Table
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors
import os
import glob
import numpy as np

dir='/beegfs/lofar/mjh/rgz/Winter/postfilter'
bits=dir.split('/')
table='post_'+bits[-2].replace('-','_') # e.g. post-Fall

os.chdir(dir)

con=mdb.connect('127.0.0.1', 'prefilter_user', 'WQ98xePI', 'prefilter', cursorclass=mdbcursors.DictCursor)

cur = con.cursor()
cur.execute('select * from %s' % table)
results=cur.fetchall()

td={}
for r in results:
    td[r['object']]=r['classification']

g=glob.glob('*.fits')
assert(len(g)==1)
t=Table.read(g[0])
rc=0
ac=0
for r in t:
    if r['Source_Name'] in td:
        if td[r['Source_Name']]==3:
            print(r['Source_Name'],'to be reviewed')
            cur.execute('update %s set classification=NULL,user=NULL where object="%s"' % (table,r['Source_Name']))
            rc+=1
    else:
        print(r['Source_Name'],'to be added')
        ac+=1

nr=0
for k in td:
    if np.sum(t['Source_Name']==k)==0:
        print(k,'to be removed from table')
        cur.execute('delete from %s where object="%s"' % (table,r['Source_Name']))
        nr+=1

print('***',ac,'new sources and',rc,'to review and',nr,'to be removed ***')

con.commit()
con.close()

