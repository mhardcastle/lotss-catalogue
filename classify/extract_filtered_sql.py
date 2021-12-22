from __future__ import print_function
from astropy.table import Table
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors

table='last'

con=mdb.connect('127.0.0.1', 'prefilter_user', 'WQ98xePI', 'prefilter', cursorclass=mdbcursors.DictCursor)

cur = con.cursor()

t=Table.read('LoTSS_DR2_v100.srl_13h.lr-full.sorted_step2_flux4.prefilter_lgz_weave_selection_last.fits')

command='select object,classification from %s' % (table)
cur.execute(command)
res=cur.fetchall()

cd={}
for r in res:
    cd[r['object']]=r['classification']

classifications=[]
for r in t:
    if r['Source_Name'] in cd:
        classifications.append(cd[r['Source_Name']])
    else:
        classifications.append(-1)

print(len(classifications), len(t))
    
t['Prefilter']=classifications

t.write('LoTSS_DR2_v100.srl_13h.lr-full.sorted_step2_flux4.prefilter_lgz_weave_selection_last_classified.fits')

tf=t[t['Prefilter']==1]
tf.write('LoTSS_DR2_v100.srl_13h.lr-full.sorted_step2_flux4.prefilter_lgz_weave_selection_last_tolgz.fits')

con.close()


            
