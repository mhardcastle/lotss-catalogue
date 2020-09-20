from astropy.table import Table
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors

table='8h'

con=mdb.connect('127.0.0.1', 'prefilter_user', 'WQ98xePI', 'prefilter', cursorclass=mdbcursors.DictCursor)

cur = con.cursor()

t=Table.read('LoTSS_DR2_v100.srl_13h.lr-full.sorted_step2_flux4.prefilter_lgz_weave_selection_2.fits')

command='select object,classification from %s' % (table)
cur.execute(command)
res=cur.fetchall()

cd={}
for r in res:
    cd[r['object']]=r['classification']

classifications=[]
for r in t:
    classifications.append(cd[r['Source_Name']])

print len(classifications), len(t)
print classifications
    
t['Prefilter']=classifications

t.write('LoTSS_DR2_v100.srl_13h.lr-full.sorted_step2_flux4.prefilter_lgz_weave_selection_1b-classified.fits')

tf=t[t['Prefilter']==1]
tf.write('LoTSS_DR2_v100.srl_13h.lr-full.sorted_step2_flux4.prefilter_lgz_weave_selection_1b-tolgz.fits')

con.close()


            
