from astropy.table import Table
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors

table='8h'

con=mdb.connect('127.0.0.1', 'prefilter_user', 'WQ98xePI', 'prefilter', cursorclass=mdbcursors.DictCursor)

cur = con.cursor()

t=Table.read('LoTSS_DR2_v100.srl_13h.lr-full.sorted_step2_flux4.prefilter_lgz_weave_selection_2.fits')
for r in t:
    command='insert into %s(object) values ("%s")' % (table,r['Source_Name'])
    cur.execute(command)

con.commit()
con.close()


            
