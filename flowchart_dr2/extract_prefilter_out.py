import sys
from astropy.table import Table, Column
import numpy as np
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors

version = 'v110'
    
h = str(sys.argv[1])
if 'h' not in h:
    h+='h'
if h not in  ['0h','13h']:
    print('unknown field code (should be 0h or 13h)',h)
    sys.exit(1)
    
    
#if h =='0h':
    #tables = ['0h']
#elif h =='13h':
    #tables =  ['13h60', '13h60b','8h','13h60fix','13h40']

tables = ['0h', '13h40', '13h60', '13h60b', '13h60fix', '8h','hetdex','last', 'post_Fall', 'post_Spring', 'post_Spring_40_45', 'post_Spring_60_65', 'post_Winter', 'mfix', 'missing','missing2']



con=mdb.connect('127.0.0.1', 'prefilter_user', 'WQ98xePI', 'prefilter', cursorclass=mdbcursors.DictCursor)

cur = con.cursor()

t=Table.read('/beegfs/lofar/wwilliams/lofar_surveys/DR2/lr/LoTSS_DR2_{version}.srl_{h}.lr-full.fits'.format(version=version,h=h))

cd={}
for table in tables:
    command='select object,classification from %s' % (table)
    cur.execute(command)
    res=cur.fetchall()
    for r in res:
        cd[r['object']]=r['classification']

has_classification = cd.keys()
classifications=[]
for r in t:
    if r['Source_Name'] in has_classification:
        classifications.append(cd[r['Source_Name']])
    else:
        classifications.append(-99)

    
t.add_column(Column(name='Prefilter', data=np.array(classifications)))
t.keep_columns(['Source_Name','Prefilter'])

t.write('/beegfs/lofar/wwilliams/lofar_surveys/DR2/LoTSS_DR2_{version}.srl_{h}.prefilter_outputs.fits'.format(version=version,h=h),overwrite=True)


con.close()


            
