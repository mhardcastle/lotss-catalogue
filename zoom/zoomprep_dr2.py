# Separate out the preparation stages (populate sql and download if
# needed) from the actual zooming for DR2.

# Note that this way we mark a list of sources for zooming by adding
# to the sql database.

from __future__ import print_function
import glob
import sys
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors
import os
from download_image_files import LofarMaps,get_legacy,get_first,get_wise

from make_catalogue import Source,make_structure,generate_table

class tzi_sql(object):
    def __init__(self,table):
        self.con=mdb.connect('127.0.0.1', 'tzi_user', 'IK34daKG', 'tzi', cursorclass=mdbcursors.DictCursor, autocommit=True)
        self.cur = self.con.cursor()
        self.execute=self.cur.execute
        self.fetchall=self.cur.fetchall
        self.table=table
        self.user=os.getenv('USER')
        sortorder=os.getenv('SORT')
        if sortorder is None:
            sortorder='ASC'
        self.sortorder=sortorder

    def get_object(self,objname,create=False):
        self.execute('select * from '+self.table+' where object=%s', (objname,))
        result=self.fetchall()
        if len(result)==0:
            if create:
                # make new record if we try to get a non-existent object
                self.execute('insert into '+self.table+' (object) values (%s)',(objname,))
                return {'object':objname}
            else:
                return None
        else:
            return result[0]

    def set_object(self,objname,record):
        for k in record:
            if k=='object':
                continue
            if record[k] is not None:
                self.execute('update '+self.table+' set '+k+'=%s where object=%s',(record[k],objname))

    def get_next(self):
        self.execute('lock table %s write' % self.table)
        self.execute('select object from %s where user is NULL and lofarfile is not NULL order by object %s limit 1' % (self.table,self.sortorder) )
        res=self.fetchall()
        if len(res)==0:
            obj=None
        else:
            obj=res[0]['object']
            self.execute('update %s set user="%s" where object="%s"' % (self.table,self.user,obj))
        self.execute('unlock tables')
        return obj

    def get_remaining(self):
        self.execute('select count(object) from %s where complete is NULL' % self.table)
        res=self.fetchall()
        remaining=int(res[0]['count(object)'])
        return remaining
            
if __name__=='__main__':

    table='Fall'
    sql=tzi_sql(table)

    s=Source.load('structure-v0.2',gaussians=False,components=False)

    os.chdir('zoom')
    sourcelist=[]
    try:
        sourcename=sys.argv[1]
    except:
        sourcename=None
    if sourcename is None:
        for sourcename in s.sd:
            if 'Deleted' in s.sd[sourcename]:
                continue
            if ( ('Zoom_prob' in s.sd[sourcename] and s.sd[sourcename]['Zoom_prob']>0.5) or
                 ('Imagemissing_prob' in s.sd[sourcename] and s.sd[sourcename]['Imagemissing_prob']>0.5)):
            
                sourcelist.append(sourcename)
    else:
        if sourcename.startswith('ILTJ'):
            sourcelist=[sourcename]
        else:
            # maybe this is a file with a list of sources
            if not os.path.isfile(sourcename):
                raise RuntimeError('Cannot parse entity on command line')
            else:
                sourcelist=[l.rstrip() for l in open(sourcename).readlines()]

    sourcelist+=s.zoomneeded
    sourcelist=list(set(sourcelist))
    print('Processing',len(sourcelist),'sources')

    # now download any images we need!
    wd=os.getcwd() # the zoom directory
    lm=LofarMaps(stay_in_imagedir=True)
    
    # now work in downloads dir
    os.chdir('downloads')
    
    for sourcename in sourcelist:
        record=sql.get_object(sourcename,create=True)
        if 'legacyfile' in record and record['legacyfile'] is not None:
            continue
        ra=s.sd[sourcename]['RA']
        dec=s.sd[sourcename]['DEC']
        record['lofarfile']=lm.find(ra,dec)
        try:
            record['legacyfile']=get_legacy(ra,dec,bands='zrg')
        except RuntimeError as e: # assume catching a 500 here
            if '500' in e.args[0]:
                print('Caught a 500, no image!')
                record['legacyfile']=None
            else:
                raise
        record['wisefile']=get_wise(ra,dec,1)
        sql.set_object(sourcename,record)

    # back to zoom dir
    os.chdir(wd)
