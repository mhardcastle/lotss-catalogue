from astropy.table import Table
from astropy.io import fits
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors
import os
import glob
import numpy as np

table='blend_filter_Spring'
gt=Table(fits.getdata('/beegfs/lofar/mjh/rgz/Spring/gaussian_lr.fits'))


#os.chdir('/beegfs/lofar/mjh/flowchart-endpoints-dr2/%s-prefilter' % table)
os.chdir('/beegfs/lofar/mjh/rgz/Spring/blend-filter')

con=mdb.connect('127.0.0.1', 'prefilter_user', 'WQ98xePI', 'prefilter', cursorclass=mdbcursors.DictCursor)
cur=con.cursor()

g=glob.glob('*.fits')
assert(len(g)==1)

t=Table.read(g[0])
for r in t:
    gaussians=gt[gt['Source_Name']==r['Source_Name']]
    all_ids=np.all(gaussians['ra']<360.0)
    command='update %s set gaussian_ids=%i where object="%s"' % (table,1 if all_ids else 0,r['Source_Name'])
    cur.execute(command)

con.commit()
con.close()


            
