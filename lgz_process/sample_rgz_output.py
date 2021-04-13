from __future__ import print_function, division
from astropy.table import Table
import sys
import numpy as np
import os
from surveys_db import SurveysDB

def run(c):
    result=os.system(c)
    if result!=0:
        raise RuntimeError('Error %i in command %s' % (result,c))

if __name__=='__main__':
    
    with SurveysDB(readonly=True) as sdb:

        sdb.cur.execute('select id,ra,decl from fields where gz_status="Complete" and dr1=0 and lgz=0')
        results=sdb.cur.fetchall()

    pos={}
    ids=[]
    for r in results:
        ids.append(r['id'])
        pos[r['id']]=(r['ra'],r['decl'])

    print('Selecting from',len(ids),'fields')

    while True:
        os.chdir('/beegfs/lofar/mjh/rgz')
        id=np.random.choice(ids)
        if os.path.isdir(id):
            print(id,'exists, skipping')
            continue
        print('Doing',id)
        ra,dec=pos[id]
        print('Position is',ra,dec)
        field_l='13h' if (ra>110 and ra<280) else '0h'
        field_h='south' if (dec<33 or field_l=='0h') else 'north'
        run('python ~/git/lotss-catalogue/lgz_process/prefilter.py filtered.csv '+id)
        os.chdir(id)
        run('python ~/git/lotss-catalogue/lgz_process/galzoo_export_deep.py %s.csv' % id)
        run('python ~/git/lotss-catalogue/lgz_process/aggregate_lofgalzoo_dr2.py %s %s' % (field_l,field_h) )
    
