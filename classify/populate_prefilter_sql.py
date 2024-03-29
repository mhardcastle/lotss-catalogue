from astropy.table import Table
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors
import os
import glob

table='blend_filter_Spring'

#os.chdir('/beegfs/lofar/mjh/flowchart-endpoints-dr2/%s-prefilter' % table)
os.chdir('/beegfs/lofar/mjh/rgz/Spring/ngp_prefilter')

con=mdb.connect('127.0.0.1', 'prefilter_user', 'WQ98xePI', 'prefilter', cursorclass=mdbcursors.DictCursor)

g=glob.glob('*.fits')
assert(len(g)==1)

cur = con.cursor()
#cur.execute('create table %s like blend_filter_Fall' % table) # 8h

t=Table.read(g[0])
for r in t:
    if os.path.isfile(r['Source_Name']+'_j.png'):
        command='insert into %s(object) values ("%s")' % (table,r['Source_Name'])
        cur.execute(command)

con.commit()
con.close()


            
