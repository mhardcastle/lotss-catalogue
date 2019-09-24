# Catalogue creation and manipulation for the deep fields

from collections import defaultdict
import os
from astropy.table import Table
import sys
import pickle
import glob
import numpy as np
from copy import deepcopy
from astropy.coordinates import SkyCoord
import astropy.units as u
import shapely
from shapely.geometry import Polygon
from shapely.ops import cascaded_union

# ellipse code taken from old process_lgz

def ellipse(x0,y0,a,b,pa,n=200):
    theta=np.linspace(0,2*np.pi,n,endpoint=False)
    st = np.sin(theta)
    ct = np.cos(theta)
    pa = np.deg2rad(pa+90)
    sa = np.sin(pa)
    ca = np.cos(pa)
    p = np.empty((n, 2))
    p[:, 0] = x0 + a * ca * ct - b * sa * st
    p[:, 1] = y0 + a * sa * ct + b * ca * st
    return Polygon(p)

class Make_Shape(object):
    '''Basic idea taken from remove_lgz_sources.py -- maybe should be merged with this one day
    but the FITS keywords are different.
    '''
    def __init__(self,clist):
        '''
        clist: a list of components that form part of the source, with RA, DEC, DC_Maj...
        '''
        ra=np.mean(clist['RA'])
        dec=np.mean(clist['DEC'])

        ellist=[]
        for r in clist:
            n_ra=r['RA']
            n_dec=r['DEC']
            x=3600*np.cos(dec*np.pi/180.0)*(ra-n_ra)
            y=3600*(n_dec-dec)
            newp=ellipse(x,y,r['DC_Maj']*3600.0+0.1,r['DC_Min']*3600.0+0.1,r['PA'])
            ellist.append(newp)
        self.cp=cascaded_union(ellist)
        self.ra=ra
        self.dec=dec
        self.h=self.cp.convex_hull
        a=np.asarray(self.h.exterior.coords)
        #for i,e in enumerate(ellist):
        #    if i==0:
        #        a=np.asarray(e.exterior.coords)
        #    else:
        #        a=np.append(a,e.exterior.coords,axis=0)
        mdist2=0
        bestcoords=None
        for r in a:
            dist2=(a[:,0]-r[0])**2.0+(a[:,1]-r[1])**2.0
            idist=np.argmax(dist2)
            mdist=dist2[idist]
            if mdist>mdist2:
                mdist2=mdist
                bestcoords=(r,a[idist])
        self.mdist2=mdist2
        self.bestcoords=bestcoords
        self.a=a

    def length(self):
        return np.sqrt(self.mdist2)

    def pa(self):
        p1,p2=self.bestcoords
        dp=p2-p1
        angle=(180*np.arctan2(dp[1],dp[0])/np.pi)-90
        if angle<-180:
            angle+=360
        if angle<0:
            angle+=180
        return angle

    def width(self):
        p1,p2=self.bestcoords
        d = np.cross(p2-p1, self.a-p1)/self.length()
        return 2*np.max(d)


def sourcename(ra,dec):
    sc=SkyCoord(ra*u.deg,dec*u.deg,frame='icrs')
    s=sc.to_string(style='hmsdms',sep='',precision=2)
    return str('ILTJ'+s).replace(' ','')[:-1]


def assemble_source(clist):
    tfluxsum=np.sum(clist['Total_flux'])
    ra=np.sum(clist['RA']*clist['Total_flux'])/tfluxsum
    dec=np.sum(clist['DEC']*clist['Total_flux'])/tfluxsum
    sname=sourcename(ra,dec)
    print 'New sourcename is',sname
    r={'Source_Name':sname}
    r['RA']=ra
    r['DEC']=dec
    r['E_RA']=np.sqrt(np.mean(clist['E_RA']**2.0))
    r['E_DEC']=np.sqrt(np.mean(clist['E_DEC']**2.0))
    r['Total_flux']=np.sum(clist['Total_flux'])
    r['E_Total_flux']=np.sqrt(np.sum(clist['E_Total_flux']**2.0))
    maxpk=np.argmax(clist['Peak_flux'])
    r['Peak_flux']=clist[maxpk]['Peak_flux']
    r['E_Peak_flux']=clist[maxpk]['E_Peak_flux']
    r['S_Code']='M'
    r['Isl_rms']=np.mean(clist['Isl_rms'])
    ms=Make_Shape(clist)
    r['LGZ_Size']=ms.length()
    r['LGZ_Width']=ms.width()
    r['LGZ_PA']=ms.pa()
    return r

class Source(object):
    def __init__(self):
        self.gd=defaultdict(dict)
        self.cd=defaultdict(dict)
        self.sd=defaultdict(dict)
        self.stage='Initialize'

    def set_stage(self,s):
        print 'Stage is',s
        self.stage=s

    def create_gaussian(self,n,r):
        # r is a table row
        for c in r.colnames:
            self.gd[n][c]=r[c]
        self.gd[n]['Created']=self.stage

    def create_component(self,n,r):
        # r is a table row
        for c in r.colnames:
            self.cd[n][c]=r[c]
        self.cd[n]['Created']=self.stage

    def create_source(self,n,r):
        # r is a table row
        try:
            keys=r.colnames
        except:
            keys=r.keys()
        for c in keys:
            self.sd[n][c]=r[c]
        self.sd[n]['Created']=self.stage

    def promote_component(self,n):
        # Promote component to source
        self.sd[n]=deepcopy(self.cd[n])
        self.sd[n]['Children']=[n]
        self.cd[n]['Parent']=n
        self.sd[n]['Created']=self.stage

    def delete_source(self,n,reason,descend=True):
        print 'Deleting source',n,'for reason',reason
        self.sd[n]['Deleted']=reason
        if descend:
            for cchild in self.sd[n]['Children']:
                if self.cd[cchild]['Parent']==n:
                    self.cd[cchild]['Deleted']=reason
                else:
                    print 'Not deleting child',cchild,'as parent mismatch'
                for gchild in self.cd[cchild]['Children']:
                    if self.gd[gchild]['Parent']==cchild:
                        self.gd[gchild]['Deleted']=reason
                    else:
                        print 'Not deleting Gaussian child',gchild,'as parent mismatch'

    def delete_gaussian(self,n,reason):
        print 'Deleting Gaussian',n,'for reason',reason
        self.gd[n]['Deleted']=reason
        
    def save(self,filename):
        '''
        Save the current state of the object.
        Parameters:
        filename -- a filename to save to
        '''
        f = file(filename, 'wb')
        pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    @staticmethod
    def load(filename):
        '''
        Load a previously saved object.
        Parameters:
        filename -- the name of a file to load
        '''
        with file(filename, 'rb') as f:
            return pickle.load(f)

def make_structure(field):
    print 'Reading data...'

    if field=='bootes':
        ct=Table.read('/beegfs/lofar/deepfields/Bootes_LR/new_fdeep_matches/Bootes_ML_RUN_fin_overlap_srl_workflow_th.fits')
        gt=Table.read('/beegfs/lofar/deepfields/Bootes_LR/new_fdeep_matches/Bootes_ML_RUN_fin_overlap_gaul_workflow_th.fits')
        preselect_dir='/beegfs/lofar/deepfields/Bootes_preselect'
        lgz_dir='/beegfs/lofar/deepfields/lgz/bootes'
        blend_dirs=['/beegfs/lofar/deepfields/Bootes_blend']
    elif field=='lockman':
        ct=Table.read('/beegfs/lofar/deepfields/Lockman_LR/updated_LR_cols/LH_ML_RUN_fin_overlap_srl_workflow_th.fits')
        gt=Table.read('/beegfs/lofar/deepfields/Lockman_LR/updated_LR_cols/LH_ML_RUN_fin_overlap_gaul_workflow_th.fits')
    elif field=='en1':
        ct=Table.read('/beegfs/lofar/deepfields/ELAIS_N1_LR/new_optcat_matches/EN1_ML_RUN_fin_overlap_srl_workflow_th.fits')
        gt=Table.read('/beegfs/lofar/deepfields/ELAIS_N1_LR/new_optcat_matches/EN1_ML_RUN_fin_overlap_gaul_workflow_th.fits')
    else:
        print 'Not in correct working directory'
        sys.exit(1)

    s=Source()

    s.set_stage('Ingest components')
    for r in ct:
        name=r['Source_Name']
        s.create_component(name,r)
        s.cd[name]['Children']=[]

    s.set_stage('Ingest Gaussians')
    for r in gt:
        name=r['Source_Name']
        s.create_gaussian(name,r)
        cname=ct[r['Source_id']]['Source_Name']
        s.gd[name]['Parent']=cname
        s.cd[cname]['Children'].append(name)

    s.set_stage('Create initial sources')
    for component in s.cd.keys():
        s.promote_component(component)

    # Apply pre-filter
    
    s.set_stage('Prefilter')
    group=[]
    source=[]
    wc=preselect_dir+'/*/workflow.txt'
    g=glob.glob(wc)
    for f in g:
        print f
        lines=open(f).readlines()
        for l in lines:
            l=l.rstrip()
            bits=l.split(',')
            group.append(int(bits[0]))
            source.append(bits[1])

    for sname,g in zip(source,group):
        s.sd[sname]['Prefilter']=g
        if g==5:
            s.delete_source(sname,'Artefact') # Artefact
        elif g==3: # Don't accept ID
            s.sd[sname]['lr_index_fin']=np.nan
            s.sd[sname]['lr_ra_fin']=np.nan
            s.sd[sname]['lr_dec_fin']=np.nan

    s.set_stage('Ingest LGZ')
    lgz_source=Table.read(lgz_dir+'/LGZ-cat.fits')
    lgz_source['Dec'].name='DEC'
    lgz_comps=Table.read(lgz_dir+'/LGZ-comps.fits')
    # create LGZ entries
    for r in lgz_source:
        sname=r['Source_Name']
        s.create_source(sname,r) # will add to what's on record for
                                 # the source, if it already exists
        s.sd[sname]['Children']=[]
        # Mark this source for LGZ assembly, which we'll do after TZI
        s.sd[sname]['LGZ_assembly_required']=True
    for r in lgz_comps:
        name=r['Comp_Name']
        sname=r['Source_Name']
        s.sd[sname]['Children'].append(name)
        s.cd[name]['Parent']=sname
        if r['Assoc']!=0:
            # remove sources that are now part of an association
            s.delete_source(name,'LGZ association',descend=False)

    # Blends
    s.set_stage('Blend processing')

    for bd in blend_dirs:
        patches=glob.glob(bd+'/ILT*.txt')
        for f in patches:
            name=f.replace(bd+'/','').replace('.txt','')
            if name not in s.sd:
                print name,'not in source list'
                sys.exit(2)
            lines=[l.rstrip() for l in open(f).readlines()]
            # parse the file
            if lines[0]=='## Flagged':
                s.delete_source(name,'Blend flagged')
            elif lines[0]=='## Unchanged':
                # source stays as it is
                s.sd[name]['Blend_prob']=0
                pass
            elif lines[0]=='## Components':
                print name,': output file to be processed!'
                components=len(s.cd[name]['Children'])
                print 'Component has',components,'Gaussians:',s.cd[name]['Children']
                # parse the component ID part
                child_ids=[]
                gaussian_names=[]
                gids=[]
                for i in range(components):
                    bits=lines[i+1].split()
                    gid=int(bits[0])
                    gname=gt[gid]['Source_Name']
                    # sanity check
                    if gname not in s.cd[name]['Children']:
                        raise RuntimeError('Gaussian %s not in child list' %gname)
                    gids.append(gid)
                    gaussian_names.append(gname)
                    child_ids.append(int(bits[1]))
                child_ids=np.array(child_ids)
                gids=np.array(gids)
                optid={}
                for l in lines[components+3:]:
                    bits=l.split()
                    id=int(bits[0])
                    ra=float(bits[1])
                    dec=float(bits[2])
                    optid[id]=(ra,dec)
                if np.all(child_ids==0):
                    print 'No unflagged components!'
                    # probably should never happen
                    s.delete_source(name,'All components removed')
                elif np.all(child_ids==1):
                    print 'Components unchanged'
                    # Using the LGZ names
                    if 1 in optid:
                        ra,dec=optid[1]
                        s.sd[name]['optRA']=ra
                        s.sd[name]['optDec']=dec
                        s.sd[name]['OptID_Name']='Altered'
                    else:
                        # opt id was removed
                        s.sd[name]['OptID_Name']="None"
                        s.sd[name]['optRA']=np.nan
                        s.sd[name]['optDec']=np.nan
                else:
                    print "It's complicated"
                    # This means that the component has been split
                    # into more than one set of Gaussians, possibly
                    # each with an optical ID.  After this processing
                    # we REMOVE the parent source and component. 
                    ss=set(child_ids)
                    for source_id in ss:
                        gaussians=[n for n,index in zip(gaussian_names, child_ids) if index==source_id]
                        if source_id==0:
                            # These Gaussians have not been included
                            # in any source. Since the whole source is
                            # not composed of artefacts, we should
                            # just be able to flag the Gaussian(s).
                            for g in gaussians:
                                s.delete_gaussian(g,'Gaussian not included in blend')
                        else:
                            if len(gaussians)==1:
                                print 'Promoting single Gaussian to component and source'
                                # Single Gaussian should be promoted to a source
                                gname=gaussians[0]
                                s.cd[gname]=deepcopy(s.gd[gname])
                                s.cd[gname]['Children']=[gname]
                                s.gd[gname]['Parent']=gname
                                s.promote_component(gname)
                                sname=gname
                            else:
                                print 'Assembling several Gaussians to component and source'
                                # Several Gaussians need to be
                                # assembled into a source, which will
                                # have a new name and other new
                                # properties. Use the Gaussian table for this
                                clist=gt[gids]
                                r=assemble_source(clist)
                                r['Assoc']=len(clist)
                                r['Assoc_Qual']=1
                                r['ID_Qual']=1
                                r['Blend_prob']=0
                                r['Created']='Deblend'
                                s.cd[sname]=deepcopy(r)
                                s.cd[sname]['Children']=gaussians
                                for g in gaussians:
                                    s.gd[g]['Parent']=sname
                                s.promote_component(sname)
                                
                            if source_id in optid:
                                ra,dec=optid[source_id]
                                s.sd[sname]['optRA']=ra
                                s.sd[sname]['optDec']=dec
                                s.sd[sname]['OptID_Name']='Altered'
                                
                    s.delete_source(name,'Removed by blend file')

                            
            else:
                raise RuntimeError('Cannot parse input file...'+lines[0])
                
    s.set_stage('Too zoomed in')

    # not written yet

    s.set_stage('LGZ post-processing')
    for name in s.sd.keys():
        if 'Deleted' not in s.sd[name].keys() and 'LGZ_assembly_required' in s.sd[name].keys():
            print 'Need to assemble source',name,'from children',s.sd[name]['Children']
            cids=[]
            for cname in s.sd[name]['Children']:
                cids.append(np.argmax(ct['Source_Name']==cname))
            clist=ct[cids]
            r=assemble_source(clist)
            r['Source_Name']=name # use the LGZ-cat sourcenames
            s.create_source(name,r)
            if s.sd[name]['Art_prob']>0.5:
                s.delete_source(name)
    

    return s

def write_table(outname,sd,columns):
    
    sources=sorted(sd.keys())
    headers=[]
    cols=[]
    for c,default in columns:
        print c,'... ',
        headers.append(c)
        column=[]
        for sn in sources:
            if 'Deleted' in sd[sn].keys():
                continue
                
            if c in sd[sn].keys():
                column.append(sd[sn][c])
            elif default is None:
                raise RuntimeError('Source %s (created by route %s) missing required key %s' % (sn,sd[sn]['Created'],c))
            else:
                column.append(default)
        cols.append(column)
    t=Table(cols,names=headers)
    t.write(outname,overwrite=True)
    

if __name__=='__main__':

    dir=os.getcwd()
    field=os.path.basename(dir)
    print 'field is',field
    s=make_structure(field)

    print 'Constructing output table'
    columns=[('Source_Name',None),('RA',None),('DEC',None),('E_RA',None),('E_DEC',None),('Total_flux',None),('E_Total_flux',None),('Peak_flux',None),('E_Peak_flux',None),('S_Code',None),('Maj',np.nan),('Min',np.nan),('PA',np.nan),('E_Max',np.nan),('E_min',np.nan),('E_PA',np.nan),('DC_Maj',np.nan),('DC_Min',np.nan),('DC_PA',np.nan),('FLAG_WORKFLOW',0),('lr_ra_fin',np.nan),('lr_dec_fin',np.nan),('optRA',np.nan),('optDec',np.nan),('LGZ_Size',np.nan),('LGZ_Width',np.nan),('LGZ_PA',np.nan),('Assoc',0),('Assoc_Qual',np.nan),('Art_prob',np.nan),('Blend_prob',np.nan),('Hostbroken_prob',np.nan),('Imagemissing_prob',np.nan),('Zoom_prob',np.nan),('Created',None)]
    write_table('test-table.fits',s.sd,columns)
            
    #s.save('test.pickle')

