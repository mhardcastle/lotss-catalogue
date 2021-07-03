from astropy.table import Table
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors
import os
import glob

dir='/beegfs/lofar/mjh/rgz/Spring-40-45/postfilter'
bits=dir.split('/')
table='post_'+bits[-2].replace('-','_') # e.g. post-Fall

os.chdir(dir)

con=mdb.connect('127.0.0.1', 'prefilter_user', 'WQ98xePI', 'prefilter', cursorclass=mdbcursors.DictCursor)

cur = con.cursor()
cur.execute('create table %s like 8h' % table)

g=glob.glob('*.fits')
assert(len(g)==1)
t=Table.read(g[0])
for r in t:
    if os.path.isfile(r['Source_Name']+'_S.png'):
        command='insert into %s(object) values ("%s")' % (table,r['Source_Name'])
        cur.execute(command)

con.commit()
con.close()


            
