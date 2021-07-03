import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors
import os

dir='/beegfs/lofar/mjh/rgz/Spring-60-65/postfilter'
bits=dir.split('/')
table='post_'+bits[-2].replace('-','_') # e.g. post-Fall
os.chdir(dir)

os.mkdir('old')

con=mdb.connect('127.0.0.1', 'prefilter_user', 'WQ98xePI', 'prefilter', cursorclass=mdbcursors.DictCursor, autocommit=True)
cur = con.cursor()
cur.execute('select * from %s where classification=3' % (table,))
results=cur.fetchall()
print('moving',len(results),'sources')
for r in results:
    os.system('mv '+r['object']+'_S.png old')

        
