# Separate out the preparation stages (populate sql and download if
# needed) from the actual zooming for DR2.

# Note that this way we mark a list of sources for zooming by adding
# to the sql database.

# Adapted version for the new_blend workflow

from __future__ import print_function
import glob
import sys
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors
import os
from download_image_files import LofarMaps,get_legacy,get_first,get_wise
from find_wise import WISE

from make_catalogue import Source,make_structure,generate_table

class prefilter_sql(object):
    def __init__(self,table):
        self.con=mdb.connect('127.0.0.1', 'prefilter_user', 'WQ98xePI', 'prefilter', cursorclass=mdbcursors.DictCursor)
        self.cur = self.con.cursor()
        self.table=table

    def get_objects(self,classification):
        command='select object from %s where classification=%i' % (self.table,classification)
        print('Running',command)
        try:
            self.cur.execute(command)
            result=self.cur.fetchall()
        except mdb.ProgrammingError as p:
            print(p)
            result=None
        if result is None:
            return result
        olist=[]
        for r in result:
            olist.append(r['object'])
        return olist
        
    def close(self):
        self.cur.close()
        self.con.close()
            
    def __exit__(self):
        close(self)
        

class blend_sql(object):
    def __init__(self,table):
        self.con=mdb.connect('192.168.2.249', 'tzi_user', 'IK34daKG', 'tzi', cursorclass=mdbcursors.DictCursor, autocommit=True)
        self.cur = self.con.cursor()
        self.execute=self.cur.execute
        self.fetchall=self.cur.fetchall
        self.table=table+'_blend'
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

    def set_object(self,objname,record,set_null=False):
        for k in record:
            if k=='object':
                continue
            if record[k] is not None or set_null:
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

    def close(self):
        self.cur.close()
        self.con.close()
            
if __name__=='__main__':

    field=os.getcwd().split('/')[-1].replace('-','_')
    table=field
    sql=blend_sql(table)

    fname=sorted(glob.glob('structure-v*-sources.pickle'))[-1].replace('-sources.pickle','')
    
    s=Source.load(fname,gaussians=False,components=False)

    sourcelist=[]
    try:
        sourcename=sys.argv[1]
    except:
        sourcename=None
    if sourcename is None:
        # no arguments, construct zoom list from source database
        for sourcename in s.sd:
            if 'Deleted' in s.sd[sourcename]:
                continue
            if 'Blend_prob' in s.sd[sourcename] and s.sd[sourcename]['Blend_prob']>0.5:
                sourcelist.append(sourcename)
        print('Starting from an auto-generated source list of length',len(sourcelist))
    else:
        if sourcename.startswith('ILTJ'):
            sourcelist=[sourcename]
        else:
            # maybe this is a file with a list of sources
            if not os.path.isfile(sourcename):
                raise RuntimeError('Cannot parse entity on command line')
            else:
                sourcelist=[l.rstrip() for l in open(sourcename).readlines()]

        print('Starting from a user-supplied source list of length',len(sourcelist))

        # check whether any of these sources have zoom files already
        for source in sourcelist:
            if source not in s.sd:
                print(source,'does not exist in source table!!')
                continue
            zoomfile=s.sd[source].get('Blend_file')
            if zoomfile and 'new_blend' in zoomfile:
                print(source,'has a blend file already:',zoomfile)
                if not os.path.isdir('old_blend'):
                    os.mkdir('old_blend')
                try:
                    os.rename(zoomfile,'old_blend/'+os.path.basename(zoomfile))
                except OSError:
                    print('... rename failed, maybe you already removed it?')
            record=sql.get_object(source,create=False)
            if record is not None:
                print(source,': clearing existing SQL record')
                sql.execute('update '+sql.table+' set complete=NULL where object="%s"' % source)
                sql.execute('update '+sql.table+' set user=NULL where object="%s"' % source)
            
    sourcelist=list(set(sourcelist))
    print('Processing',len(sourcelist),'sources')

    os.chdir('new_blend')
    # now download any images we need!
    wd=os.getcwd() # the zoom directory
    lm=LofarMaps(stay_in_imagedir=True)
    w=WISE()
    
    # now work in downloads dir
    os.chdir('downloads')
    
    for sourcename in sourcelist:
        print(sourcename)
        record=sql.get_object(sourcename,create=True)
        if 'legacyfile' in record and record['legacyfile'] is not None:
            continue
        if sourcename not in s.sd:
            print('Skipping source',sourcename,'with no record')
            continue
        if 'RA' not in s.sd[sourcename]:
            print('Skipping source',sourcename,'with no RA in record: record follows')
            print(s.sd[sourcename])
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
        wisename=w.find_pos(ra,dec)
        if wisename is None:
            wisename=get_wise(ra,dec,1)
        record['wisefile']=wisename
        sql.set_object(sourcename,record)

    # back to zoom dir
    os.chdir(wd)
