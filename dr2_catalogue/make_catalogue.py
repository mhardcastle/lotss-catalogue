# Catalogue creation and manipulation for DR2 RGZ

from __future__ import print_function
from collections import defaultdict
import os
from astropy.table import Table
import sys
import cPickle as pickle
import glob
import numpy as np
from copy import deepcopy
from astropy.coordinates import SkyCoord
import astropy.units as u
import shapely
from shapely.geometry import Polygon
from shapely.ops import cascaded_union
from separation import separation
from tqdm import tqdm
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors

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
            newp=ellipse(x,y,r['DC_Maj']+0.1,r['DC_Min']+0.1,r['PA'])
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
    #print('New sourcename is',sname)
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
    r['S_Code']='Z'
    r['Mosaic_ID']=clist[0]['Mosaic_ID']
    r['Isl_rms']=np.mean(clist['Isl_rms'])
    ms=Make_Shape(clist)
    r['LGZ_Size']=ms.length()
    #if r['LGZ_Size']>1000:
    #    print(clist['Source_Name'])
    #    print('Unreasonable source size detected')
    r['LGZ_Width']=ms.width()
    r['LGZ_PA']=ms.pa()
    for k in ['Maj','Min','PA','E_Maj','E_Min','E_PA','DC_Maj','DC_Min','DC_PA']:
        r[k]=np.nan
    return r

def parsefile(sourcename,ss,dir=''):
    lines=[l.rstrip() for l in open(dir+sourcename+'.txt').readlines()]
    ss.sd[sourcename]['Zoom_prob']=0
    ss.sd[sourcename]['Imagemissing_prob']=0
    ss.sd[sourcename]['Hostbroken_prob']=0
    if 'Children' not in ss.sd[sourcename]:
        # because this is a zoom file for a source that doesn't exist
        ss.sd[sourcename]['Children']=[]
    if 'Deleted' in ss.sd[sourcename]:
        del(ss.sd[sourcename]['Deleted'])
    lc=0
    ra=None
    dec=None
    while lc<len(lines):
        l=lines[lc]
        if l=="":
            lc+=1
            continue
        if l[0]=='#':
            bits=l.split()
            if bits[1]=="Components":
                comps=[]
                lc+=1
                while lc<len(lines) and lines[lc]!="":
                    if ss.iscomp(lines[lc]):
                        comps.append(lines[lc])
                    lc+=1
                ss.set_components(sourcename,comps,flag_removals=True)
            elif bits[1]=="OptID":
                lc+=1
                l=lines[lc]
                ra,dec=[float(b) for b in l.split()]
                if ra<=360 and dec<=90:
                    ss.set_opt(sourcename,ra,dec)
                else:
                    ss.set_opt(sourcename,np.nan,np.nan)
                lc+=1
            elif bits[1]=="Size":
                ss.set_size(sourcename,float(lines[lc+1]))
                lc+=1
            elif bits[1]=="Deleted":
                ss.delete_source(sourcename,'Zoom marked artefact')
                lc+=1
                while lc<len(lines) and lines[lc]!="":
                    csource=lines[lc]
                    if csource in ss.sd and 'Deleted' not in ss.sd[csource]:
                        ss.delete_source(csource,'Zoom marked artefact')
                    # this shouldn't be necessary, but just in case
                    if ss.iscomp(csource) and 'Deleted' not in ss.cd[csource]:
                        ss.delete_component(csource,'Zoom marked artefact')
                    lc+=1
            elif bits[1]=="Blend":
                ss.sd[sourcename]['Blend_prob']=1.0
        lc+=1
    if ra is None:
        ss.set_opt(sourcename,np.nan,np.nan)
        ss.sd[sourcename]['NoID']=11
        
class Source(object):
    def __init__(self,logfile=None):
        self.gd=defaultdict(dict)
        self.cd=defaultdict(dict)
        self.sd=defaultdict(dict)
        self.zoomneeded=[]
        self.zoomreasons=[]
        self.stage='Initialize'
        self.logfile=logfile

    def log(self,*kwargs):
        print(*kwargs,file=self.logfile)
        
    def set_stage(self,s):
        print('Stage is',s)
        self.log('=========== Stage is',s,' ============')
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
        if 'Created' not in self.sd[n]:
            self.sd[n]['Created']=self.stage

    def promote_component(self,n):
        # Promote component to source
        self.sd[n]=deepcopy(self.cd[n])
        self.sd[n]['Children']=[n]
        self.cd[n]['Parent']=n
        self.sd[n]['Created']=self.stage

    def delete_source(self,n,reason,descend=True):
        self.log('Deleting source',n,'for reason',reason)
        self.sd[n]['Deleted']=reason
        if descend:
            for cchild in self.sd[n]['Children']:
                if self.cd[cchild]['Parent']==n:
                    self.cd[cchild]['Deleted']=reason
                else:
                    self.log('Not deleting child',cchild,'as parent mismatch')

                for gchild in self.cd[cchild]['Children']:
                    if self.gd[gchild]['Parent']==cchild:
                        self.gd[gchild]['Deleted']=reason
                    else:
                        self.log('Not deleting Gaussian child',gchild,'as parent mismatch')

    def delete_component(self,n,reason,descend=True):
        self.log('Deleting component',n,'for reason',reason)
        self.cd[n]['Deleted']=reason
        for gchild in self.cd[n]['Children']:
            if self.gd[gchild]['Parent']==n:
                self.gd[gchild]['Deleted']=reason
            else:
                self.log('Not deleting Gaussian child',gchild,'as parent mismatch')
                       
    def delete_gaussian(self,n,reason):
        self.log('Deleting Gaussian',n,'for reason',reason)
        self.gd[n]['Deleted']=reason

    # methods for the zoom code

    def get_comps(self,sourcename):
        components=[]
        if 'Children' in self.sd[sourcename]:
            for c in self.sd[sourcename]['Children']:
                if c in self.cd and 'Deleted' not in self.cd[c]:
                    components.append(c)
        else:
            self.log('Source',sourcename,'has no children -- this should not happen!')
        return components

    def get_ncomps(self,sourcename):
        ncomp=[]
        for k in self.sd:
            if k!=sourcename and 'Deleted' not in self.sd[k]:
                if 'Children' not in self.sd[k]:
                    self.log('Children missing!',k,self.sd[k])
                else:
                    ncomp+=self.sd[k]['Children']
        return ncomp

    def addzoom(self,sourcename,reason):
        self.log('Adding',sourcename,'to zoom list for reason:',reason)
        self.zoomneeded.append(sourcename)
        self.zoomreasons.append(reason)
    
    def set_components(self,sourcename,componentlist,flag_removals=False):
        self.log('setting components of',sourcename,'to',componentlist)
        orphans=list(set(self.sd[sourcename]['Children'])-set(componentlist))
        self.sd[sourcename]['Children']=componentlist
        for comp in componentlist:
            # do these components have other parents?
            if 'Parent' in self.cd[comp] and self.cd[comp]['Parent']!='':
                oldparent=self.cd[comp]['Parent']
                if oldparent!=sourcename: # parent has changed
                    if flag_removals and oldparent!=comp:
                        # must be a previous LGZ source
                        self.addzoom(oldparent,'set_components parent changed')
                    if 'Children' not in self.sd[oldparent]:
                        self.log('Old parent source of component %s (%s) has no children!' % (comp, oldparent))
                        self.log(self.cd[comp])
                        self.delete_source(oldparent,'Reallocated (orphan)')
                    else:
                        self.sd[oldparent]['Children'].remove(comp)
                        if self.sd[oldparent]['Children']==[]:
                            # remove broken flag if source now has no children
                            if 'Broken' in self.sd[oldparent]:
                                del(self.sd[oldparent]['Broken'])
                            self.delete_source(oldparent,'Reallocated (no children)')
                        else:
                            # flag a composite source that has had
                            # children removed. These should be
                            # examined
                            self.sd[oldparent]['Broken']='Composite source has had children removed'

            self.cd[comp]['Parent']=sourcename
            #over-ride component deletion if a zoom file specifies it
            if 'Deleted' in self.cd[comp]:
                del(self.cd[comp]['Deleted'])
            '''
            # skip check of all sources
            for checksource in self.sd:
                if 'Children' not in self.sd[checksource]:
                    self.log('Children field missing!',checksource,self.sd[checksource])
                    self.sd[checksource]['Children']=[]
                if checksource!=sourcename and comp in self.sd[checksource]['Children']:
                   self.sd[checksource]['Children'].remove(comp)
                   if self.sd[checksource]['Children']==[]:
                       self.delete_source(checksource,'Reallocated in zoom')
            '''
        for comp in orphans:
            self.log('*** Creating orphan %s' % comp)
            self.cd[comp]['Parent']=''
            
    def set_opt(self,sourcename,ra,dec):
        self.sd[sourcename]['optRA']=ra
        self.sd[sourcename]['optDec']=dec
        self.sd[sourcename]['lr_ra_fin']=np.nan
        self.sd[sourcename]['lr_dec_fin']=np.nan

    def set_size(self,sourcename,size):
        self.sd[sourcename]['Manual_Size']=size
    
    def optid(self,sourcename):
        if 'optRA' in self.sd[sourcename]:
            return self.sd[sourcename]['optRA'],self.sd[sourcename]['optDec']
        elif 'lr_ra_fin' in self.sd[sourcename]:
            return self.sd[sourcename]['lr_ra_fin'],self.sd[sourcename]['lr_dec_fin']
        else:
            return np.nan,np.nan
            
    def iscomp(self,component):
        return component in self.cd
        
    def save(self,filename):
        '''
        Save the current state of the object.
        Parameters:
        filename -- a base filename to save to
        '''
        with open(filename+'-sources.pickle', 'wb') as f:
            pickle.dump(self.sd, f, pickle.HIGHEST_PROTOCOL)
        with open(filename+'-components.pickle', 'wb') as f:
            pickle.dump(self.cd, f, pickle.HIGHEST_PROTOCOL)
        with open(filename+'-gaussians.pickle', 'wb') as f:
            pickle.dump(self.gd, f, pickle.HIGHEST_PROTOCOL)
        with open(filename+'-zoomneeded.pickle', 'wb') as f:
            pickle.dump(self.zoomneeded, f, pickle.HIGHEST_PROTOCOL)
        with open(filename+'-zoomreasons.pickle', 'wb') as f:
            pickle.dump(self.zoomreasons, f, pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def load(filename,components=True,gaussians=True):
        '''
        Load a previously saved object.
        Parameters:
        filename -- the base name of files to load
        '''
        s=Source() # an empty instance
        s.stage='Load'
        print('Loading sources...')
        with open(filename+'-sources.pickle', 'rb') as f:
            s.sd=pickle.load(f)
        if components:
            print('Loading components...')
            with open(filename+'-components.pickle', 'rb') as f:
                s.cd=pickle.load(f)
        if gaussians:
            print('Loading Gaussians...')
            with open(filename+'-gaussians.pickle', 'rb') as f:
                s.gd=pickle.load(f)
        
        with open(filename+'-zoomneeded.pickle', 'rb') as f:
            s.zoomneeded=pickle.load(f)

        with open(filename+'-zoomreasons.pickle', 'rb') as f:
            s.zoomreasons=pickle.load(f)
        return s

def warn_or_die(warn,s):
    if warn:
        print('*** WARNING: '+s)
    else:
        raise RuntimeError(s)
        
def make_structure(field,warn=False,version=None):
    print('Reading data...')

    lgz_dir='.'
    ct=Table.read('source_lr.fits')
    ct['ra'].name='lr_ra_fin'
    ct['dec'].name='lr_dec_fin'
    ct['lr_ra_fin']=np.where(ct['lr_ra_fin']>360,np.nan,ct['lr_ra_fin'])
    ct['lr_dec_fin']=np.where(ct['lr_dec_fin']>90,np.nan,ct['lr_dec_fin'])
    ct['lr'].name='lr_fin'
    gt_outname=lgz_dir+'/edited_gaussians.fits'
    if os.path.isfile(gt_outname):
        gt=Table.read(gt_outname)
    else:
        gt=Table.read('gaussian_lr.fits')
        print('Fixing up Gaussian table')
        gt['Source_Name'].name='Parent_Name'
        sc=SkyCoord(gt['RA'],gt['DEC'],frame='icrs')
        strings=sc.to_string(style='hmsdms',sep='',precision=2)
        ilt=[]
        for s in strings:
            ilt.append(str('ILTJ'+s).replace(' ','')[:-1])
        gt['Source_Name']=ilt

        gt['Gaus_id']=list(range(len(gt)))
        print('Writing amended file')
        gt.write(gt_outname)
        
    preselect_dir=None
    blend_dirs=['blend']
    noid_files=[]
    ridgeline='allhosts.fits'
    if version is None:
        logfilename='logfile.txt'
    else:
        logfilename='logfile-'+version+'.txt'
    
    s=Source(logfile=open(logfilename,'w'))
    
    s.set_stage('Ingest components')
    for r in tqdm(ct):
        name=r['Source_Name']
        s.create_component(name,r)
        s.cd[name]['Children']=[]

    s.set_stage('Ingest Gaussians')
    for r in tqdm(gt):
        name=r['Source_Name']
        cname=r['Parent_Name']
        if cname not in s.cd:
            continue # drop orphan Gaussian
        s.create_gaussian(name,r)
        s.gd[name]['Parent']=cname
        s.cd[cname]['Children'].append(name)

    s.set_stage('Create initial sources')
    artefacts=[]
    for component in tqdm(s.cd):
        s.promote_component(component)
        if s.cd[component]['Prefilter']==4:
            s.addzoom(component,'Prefilter zoom required')
        if s.cd[component]['Prefilter']==5:
            artefacts.append(component)

    print('Deleting',len(artefacts),'prefilter artefacts')
    for c in tqdm(artefacts):
        s.delete_source(c,'Artefact')
            
    s.set_stage('Ingest LGZ')
    lgz_source=Table.read(lgz_dir+'/LGZ-cat.fits')
    lgz_source['Dec'].name='DEC'
    lgz_comps=Table.read(lgz_dir+'/LGZ-comps.fits')
    # create LGZ entries
    for r in tqdm(lgz_source):
        sname=r['Source_Name']
        s.create_source(sname,r) # will add to what's on record for
                                 # the source, if it already exists
        if s.sd[sname]['Compoverlap']==1:
            s.addzoom(sname,'Compoverlap set')
        if 'lr_ra_fin' in s.sd[sname]:
            s.sd[sname]['old_ra']=s.sd[sname]['lr_ra_fin']
            s.sd[sname]['old_dec']=s.sd[sname]['lr_dec_fin']
        s.sd[sname]['lr_ra_fin']=np.nan
        s.sd[sname]['lr_dec_fin']=np.nan
        s.sd[sname]['Children']=[]
        # Mark this source for LGZ assembly, which we'll do after TZI
        s.sd[sname]['LGZ_assembly_required']=True
        # make a component table
        source_comps=lgz_comps[lgz_comps['Source_Name']==sname]
        s.set_components(sname,list(source_comps['Comp_Name']),flag_removals=True)
        '''
    for r in lgz_comps:
        name=r['Comp_Name']
        sname=r['Source_Name']
        if name not in s.cd:
            raise RuntimeError('Component %s from source %s in LGZ does not exist in table' % (name,sname))
        s.sd[sname]['Children'].append(name)
        s.cd[name]['Parent']=sname
        if r['Assoc']!=0:
            # remove sources that are now part of an association
            s.delete_source(name,'LGZ association',descend=False)
        '''

    # Blends
    s.set_stage('Blend processing')

    for bd in blend_dirs:
        patches=glob.glob(bd+'/ILT*.txt')
        #if len(patches)==0:
        #    raise RuntimeError('No patches found in directory %s, did you get the pathname wrong?' % bd)
        for f in tqdm(patches):
            name=f.replace(bd+'/','').replace('.txt','')
            s.log('Blend file for',name,'is',f)
            if name not in s.sd:
                raise RuntimeError(name+' not in source list')
            # whatever happens to the source, these issues are dealt with...
            s.sd[name]['Blend_prob']=0
            s.sd[name]['Hostbroken_prob']=0
            s.sd[name]['Imagemissing_prob']=0
            s.sd[name]['Blend_file']=f
            lines=[l.rstrip() for l in open(f).readlines()]
            # parse the file
            if lines[0]=='## Flagged':
                s.sd[name]['Zoom_prob']=1.0 # send to TZI
            elif lines[0]=='## Unchanged':
                # source stays as it is, BUT lr must be accepted if present
                if not np.isnan(s.sd[name]['lr_ra_fin']):
                    s.log('Accepting LR position for this source')
                    s.sd[name]['optRA']=s.sd[name]['lr_ra_fin']
                    s.sd[name]['optDec']=s.sd[name]['lr_dec_fin']
                elif 'old_ra' in s.sd[name]:
                    s.log('Accepting pre-LGZ LR position for this source')
                    s.sd[name]['optRA']=s.sd[name]['old_ra']
                    s.sd[name]['optDec']=s.sd[name]['old_dec']
                    
            elif lines[0]=='## Components':
                s.log(name,': output file to be processed!')
                s.sd[name]['lr_ra_fin']=np.nan
                s.sd[name]['lr_dec_fin']=np.nan
                components=len(s.cd[name]['Children'])
                s.log('Component has',components,'Gaussians:',s.cd[name]['Children'])
                # parse the component ID part
                child_ids=[]
                gaussian_names=[]
                gids=[]
                for i in range(components):
                    bits=lines[i+1].split()
                    gid=int(bits[0])
                    gidi=np.argmax(gt['Gaus_id']==gid)
                    gname=gt[gidi]['Source_Name']
                    # sanity check
                    if gname not in s.cd[name]['Children']:
                        print(s.cd[name])
                        raise RuntimeError('Gaussian %s not in child list' %gname)
                    gids.append(gidi)
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
                    s.log('No unflagged components!')
                    # probably should never happen
                    s.delete_source(name,'All components removed')
                elif np.all(child_ids==1):
                    s.log('Components unchanged')
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
                        s.sd[name]['noID']=11
                else:
                    s.log("It's complicated")
                    # This means that the component has been split
                    # into more than one set of Gaussians, possibly
                    # each with an optical ID.
                    # we delete the parent source now -- if it needs to be recreated it will be.
                    s.delete_source(name,'Removed by blend file',descend=False)
                    
                    ss=set(child_ids)
                    for source_id in ss:
                        gaussians=[n for n,index in zip(gaussian_names, child_ids) if index==source_id]
                        these_gids=gids[child_ids==source_id]
                        if source_id==0:
                            # These Gaussians have not been included
                            # in any source. Since the whole source is
                            # not composed of artefacts, we should
                            # just be able to flag the Gaussian(s).
                            for g in gaussians:
                                s.delete_gaussian(g,'Gaussian not included in blend')
                        else:
                            if len(gaussians)==1:
                                s.log('Promoting single Gaussian to component and source')
                                # Single Gaussian should be promoted to a source
                                gname=gaussians[0]
                                parent=s.gd[gname]['Parent']
                                if 'Promoted' not in s.cd[parent]['Created']:
                                    s.log(gname,'has parent',parent,' -- marking deleted')
                                    s.cd[parent]['Deleted']='Removed by blend file (single)'
                                else:
                                    s.log('Not deleting parent',parent)
                                s.cd[gname]=deepcopy(s.gd[gname])
                                s.cd[gname]['Created']='Promoted from single Gaussian'
                                s.cd[gname]['Children']=[gname]
                                s.gd[gname]['Parent']=gname
                                s.promote_component(gname)
                                sname=gname
                            else:
                                s.log('Assembling several Gaussians to component and source')
                                # Several Gaussians need to be
                                # assembled into a source, which will
                                # have a new name and other new
                                # properties. Use the Gaussian table for this
                                # ** NEW ** promote all Gaussians to components
                                clist=gt[these_gids]
                                r=assemble_source(clist)
                                sname=r['Source_Name']
                                while sname in s.sd:
                                    # set zoom prob to 1, because
                                    # there must be a source matching
                                    # this new one in position
                                    r['Zoom_prob']=1 
                                    if sname[-1].isdigit():
                                        sname=sname+'a'
                                    else:
                                        sname=sname[:-1]+chr(ord(sname[-1])+1)
                                    r['Source_Name']=sname
                                    s.log('Trying new source name',sname)
                                if len(clist)==1:
                                    r['Assoc']=0
                                else:
                                    r['Assoc']=len(clist)
                                r['Assoc_Qual']=1
                                r['ID_Qual']=1
                                r['Blend_prob']=0
                                r['Created']='Deblend'
                                r['Blend_file']=f
                                r['Children']=gaussians
                                s.log('Creating source',sname,'with children',gaussians)
                                #print(clist['Source_Name'])
                                #if sname in s.sd:
                                #    print('Source exists! Previous source was created by',s.sd[sname]['Created'],'and has children',s.sd[sname]['Children'])
                                #    if s.sd[sname]['Created']=='Deblend':
                                #        print('Previous blend file was',s.sd[sname]['Blend_file'],': this file is',f)
                                #    raise RuntimeError('Duplicate source name',sname)
                                s.sd[sname]=r
                                for g in gaussians:
                                    # logic here and above deals with the case where a component created elsewhere in the loop has the same name as a parent component of a Gaussian that would normally be deleted. As the old parent has already been overwritten the deletion is not necessary.
                                    parent=s.gd[g]['Parent']
                                    if 'Promoted' not in s.cd[parent]['Created']:
                                        s.log(g,'has parent',parent,' -- marking deleted')
                                        s.cd[parent]['Deleted']='Removed by blend file (multiple)'
                                    else:
                                        s.log('Not deleting parent',parent)
                                    s.cd[g]=deepcopy(s.gd[g])
                                    s.cd[g]['Created']='Promoted from Gaussian'
                                    s.cd[g]['Deblended_from']=name
                                    s.gd[g]['Parent']=g
                                    s.cd[g]['Children']=[g]
                                    s.cd[g]['Parent']=sname

                            s.sd[sname]['lr_ra_fin']=np.nan
                            s.sd[sname]['lr_dec_fin']=np.nan
                            if source_id in optid:
                                ra,dec=optid[source_id]
                                s.sd[sname]['optRA']=ra
                                s.sd[sname]['optDec']=dec
                                s.sd[sname]['OptID_Name']='Altered'
                                s.sd[sname]['NoID']=0
                            else:
                                s.sd[sname]['NoID']=11
                                s.sd[sname]['OptID_Name']="None"
                                s.sd[sname]['optRA']=np.nan
                                s.sd[sname]['optDec']=np.nan
                            
            else:
                raise RuntimeError('Cannot parse input file...'+lines[0])

    s.set_stage('Too zoomed in')

    g=sorted(glob.glob(lgz_dir+'/zoom/ILTJ*.txt'),key=os.path.getmtime)
    for f in tqdm(g):
        s.log('Zoomfile',f)
        source=f.replace('.txt','').replace(lgz_dir+'/zoom/','')
        parsefile(source,s,dir=lgz_dir+'/zoom/')
        s.sd[source]['Created']='Too zoomed in'
        s.sd[source]['Zoomfile']=f
        s.sd[source]['LGZ_assembly_required']=True

    # component table won't now change, so generate it so it can be
    # passed to assemble_source
    s.set_stage('Building new component table')
    columns=[('Source_Name',None),('RA',None),('DEC',None),('E_RA',None),('E_DEC',None),('Total_flux',None),('E_Total_flux',None),('Peak_flux',None),('E_Peak_flux',None),('S_Code',None),('Mosaic_ID',None),('Maj',np.nan),('Min',np.nan),('PA',np.nan),('E_Maj',np.nan),('E_Min',np.nan),('E_PA',np.nan),('DC_Maj',np.nan),('DC_Min',np.nan),('DC_PA',np.nan),('Isl_rms',np.nan),('Created',None),('Parent',None)]
    new_ct=generate_table(s.cd,columns)
    # create a lookup dictionary so we don't have to inefficiently search for source name many times
    print('Build the lookup dictionary for this table')
    idd={}
    for i,r in tqdm(enumerate(new_ct),total=len(new_ct)):
        idd[r['Source_Name']]=i
    
    for assembly_stage in [1,2]:
        s.set_stage('Source assembly stage %i' % assembly_stage)
        sources=s.sd.keys() # copy because we rename as we go
        for name in tqdm(sources):
            if 'Deleted' not in s.sd[name] and 'LGZ_assembly_required' in s.sd[name]:
                if len(s.sd[name]['Children'])==1 and s.sd[name]['Children'][0]==name:
                    #print('Source',name,'marked as assembly required but has only one component')
                    del s.sd[name]['LGZ_assembly_required']
                    continue # source is single component
                s.log('Need to assemble source',name,'from children',s.sd[name]['Children'])
                cids=[]
                error=False
                for cname in s.sd[name]['Children']:
                    try:
                        cid=idd[cname]
                    except KeyError:
                        cid=None
                    if cid is None:
                        s.log('Source created by',s.sd[name]['Created'])
                        error=True
                        s.log('Child %s does not exist!' % cname)
                    else:
                        cids.append(cid)
                if len(cids)==0:
                    s.delete_source(name,'All components removed')
                    continue
                if error:
                    s.log('Source is partial, needs zoom file fix')
                    s.addzoom(name,'Partial source in LGZ post-processing')
                clist=new_ct[cids]
                r=assemble_source(clist)
                if 'Manual_Size' in s.sd[name]:
                    s.log('Adding manual size measurement')
                    r['LGZ_Size']=s.sd[name]['Manual_Size']
                sname=r['Source_Name']
                if sname!=name:
                    s.log('Renaming old source',name,'created by',s.sd[name]['Created'],'to',sname)
                    r['Renamed_from']=name
                    s.delete_source(name,'Renamed',descend=False)


                for key in s.sd[name]:
                    if key not in r:
                        r[key]=s.sd[name][key] # includes copying children
                for comp in r['Children']:
                    s.cd[comp]['Parent']=sname

                '''
                # this part should only be necessary if Parent keys are screwed up
                if name!=sname:
                    components=list(s.cd)
                    for comp in components:
                        if s.cd[comp]['Parent']==name:
                            s.delete_component(comp,'Orphaned')
                '''

                r['Assoc']=len(r['Children'])
                s.create_source(sname,r)
                if 'Deleted' in s.sd[sname]:
                    del(s.sd[sname]['Deleted'])
                if 'Art_prob' in s.sd[sname] and s.sd[sname]['Art_prob']>0.5:
                    s.delete_source(sname,'LGZ artefact')
                del s.sd[sname]['LGZ_assembly_required']
        if assembly_stage==1:
            s.set_stage('Duplicates')
            # Directly read duplicate information from the database
            con=mdb.connect('192.168.2.249', 'tzi_user', 'IK34daKG', 'duplicates', cursorclass=mdbcursors.DictCursor, autocommit=True)
            cur = con.cursor()
            table=field.replace('-','_')
            try:
                cur.execute('select * from %s where classification is not NULL' % table)
                results=cur.fetchall()
            except mdb.ProgrammingError:
                results=None
            con.close()
            if results is not None:
                for r in tqdm(results):
                    # options are
                    # 1. merge sources
                    # 2. drop source 2
                    # 3. drop source 1
                    # 4. pass to TZI
                    # 5. drop optical ID (for both) and merge
                    # 6. drop optical ID (for both)  and don't merge
                    if r['source1'] not in s.sd:
                        s.log('%s does not exist in dedupe!' % r['source1'])
                        continue
                    if r['source2'] not in s.sd:
                        s.log('%s does not exist in dedupe!' % r['source2'])
                        continue
                    if r['classification']==1 or r['classification']==5:
                        # merger
                        # we take all components from source 2 and add them to source 1,
                        # then delete source 2 and mark source 1 as needing assembly.
                        s.set_components(r['source1'],s.sd[r['source1']]['Children']+s.sd[r['source2']]['Children'])
                        s.sd[r['source1']]['Created']='Deduplicate'
                        s.sd[r['source1']]['LGZ_assembly_required']=True
                        if r['classification']==5:
                            s.set_opt(r['source1'],np.nan,np.nan)
                    elif r['classification']==2:
                        s.delete_source(r['source2'],'Deleted in dedupe')
                    elif r['classification']==3:
                        s.delete_source(r['source1'],'Deleted in dedupe')
                    elif r['classification']==4:
                        s.addzoom(r['source1'],'Flagged in dedupe')
                        s.addzoom(r['source2'],'Flagged in dedupe')
                    elif r['classification']==6:
                        s.set_opt(r['source1'],np.nan,np.nan)
                        s.set_opt(r['source2'],np.nan,np.nan)
                    else:
                        raise RuntimeError('unexpected classification')

            else:
                print('No duplicates table, hope this is ok')
    # finally sort out optical positions
    s.set_stage('Sorting optical positions')
    for source in tqdm(s.sd):
        if 'Deleted' in s.sd[source]:
            continue
        if 'optRA' in s.sd[source]:
            # RH 'or' because we may have a TZI source with no opt ID
            s.sd[source]['Position_from']='Visual inspection'
        elif 'lr_ra_fin' in s.sd[source]:
            s.sd[source]['Position_from']='LR'
            s.sd[source]['optRA']=s.sd[source]['lr_ra_fin']
            s.sd[source]['optDec']=s.sd[source]['lr_dec_fin']
        else:
            s.sd[source]['Position_from']='None'

    if noid_files is not None:
        s.set_stage('NoID')
        for noid_file in noid_files:
            lines=open(noid_file).readlines()
            group=[]
            source=[]
            for l in lines:
                l=l.rstrip()
                bits=l.split(',')
                group.append(int(bits[0]))
                source.append(bits[1])
            for sname,g in zip(source,group):
                if sname not in s.sd:
                    s.log('Source',name,'already deleted, skipping')
                    continue
                #if 'optRA' in s.sd[sname] and not np.isnan(s.sd[sname]['optRA']):
                #    print 'Source',sname,'in noid list but has id, skipping'
                #    continue 
                s.sd[sname]['NoID']=g
                if g==6:
                    s.delete_source(sname,'Artefact') # Artefact
                elif g==8:
                    if 'Zoomfile' not in s.sd[sname]:
                        if 'Renamed_from' in s.sd[sname]:
                            name=s.sd[sname]['Renamed_from']
                        else:
                            name=sname
                        s.log('Adding',name,'to zoom list')
                        s.addzoom(name,'NoID file required zoom')
                    else:
                        s.sd[sname]['NoID']=3
                elif g==7 or g==9:
                    if 'Zoomfile' not in s.sd[sname]:
                        if 'LGZ' in s.sd[sname]['Created']:
                            if 'Renamed_from' in s.sd[sname]:
                                name=s.sd[sname]['Renamed_from']
                            else:
                                name=sname
                            s.log('Adding',name,'to zoom list')
                            s.addzoom(name,'NoID file required zoom (2)')
                    else:
                        s.sd[sname]['NoID']=3

    s.set_stage('Ridge line ingest')
    if ridgeline is not None and os.path.isfile(ridgeline):
        tr=Table.read(ridgeline)
        for r in tqdm(tr):
            if r['LRMagBoth']<1.0:
                continue
            name=r['Source_Name']
            if name not in s.sd:
                s.log('Warning: source %s in ridgeline file does not exist' % name)
            else:
                if 'optRA' not in s.sd[name] or np.isnan(s.sd[name]['optRA']) or 'Position_from' not in s.sd[name] or s.sd[name]['Position_from']=='LR':
                    s.sd[name]['optRA']=r['optRA_RLC']
                    s.sd[name]['optDec']=r['optDEC_RLC']
                    s.sd[name]['Position_from']='Ridge line code'
                    s.sd[name]['lr_fin']=r['LRMagBoth']

    s.set_stage('Opt ID overlap check')

    names=[]
    optras=[]
    optdecs=[]
    for name in s.sd:
        if 'Deleted' in s.sd[name]:
            continue
        if 'optRA' in s.sd[name]:
            names.append(name)
            optras.append(s.sd[name]['optRA'])
            optdecs.append(s.sd[name]['optDec'])
    names=np.array(names)
    optras=np.array(optras)
    optdecs=np.array(optdecs)
    with open('optical.pickle','w') as pf:
        pickle.dump((names,optras,optdecs),pf)
    path=os.environ['LGZPATH']
    command=path+'/dr2_catalogue/process_overlap.py'
    print('Running',command)
    retval=os.system(command)
    if retval!=0:
        raise RuntimeError('Failed to run overlap code, return value is %i' % retval)
    with open('badlist.pickle') as pf:
        bad=pickle.load(pf)
    
    badlist=[]
    badzoom=[]
    print('Processing bad sources')
    for sname in tqdm(bad):
        s.log('Checking duplicate ID source',sname,s.sd[sname]['Created'])
        if 'Zoomfile' not in s.sd[sname]:
            s.log('Adding',sname,'to zoom list')
            badlist.append(sname)
        else:
            s.log('*** problem -- source',sname,'with zoom file',s.sd[sname]['Zoomfile'],'has duplicate ID')
            badzoom.append((sname,s.sd[sname]['Zoomfile']))
                    

    badlist=list(set(badlist))
    print('Adding',len(badlist),'unique duplicate non-zoom sources to zoom list')
    for name in badlist:
        s.addzoom(name,'Duplicate ID non-zoom source')

    print('Writing list of bad zoom sources to text file for you -- check zoom files manually')
    with open('badzoom.txt','w') as f:
        for source,zoomfile in badzoom:
            f.write('%s %s\n' % (source,zoomfile))

    print('Now writing you a list of sources you have broken with excessive zoom files, some may be the same ones!')
    with open('broken.txt','w') as f:
        for source in s.sd:
            if 'Broken' in s.sd[source]:
                if 'Zoomfile' in s.sd[source]:
                    zoomfile=s.sd[source]['Zoomfile']
                else:
                    zoomfile=""
                f.write('%s %s %s\n' % (source,s.sd[source]['Created'],zoomfile))
    
    s.set_stage('Check rescue file')
    if os.path.isfile('rescue.txt'):
        lines=open('rescue.txt').readlines()
        for l in tqdm(lines):
            bits=l.rstrip().split()
            dist=float(bits[2])
            if dist>3:
                cname=bits[0]
                if cname not in s.cd:
                    raise RuntimeError('Failed to find rescued component')
                sname=s.cd[cname]['Parent']
                if sname not in s.sd:
                    s.log('Rescue source',sname,'not in source list')
                else:
                    if 'Zoomfile' not in s.sd[sname]:
                        name=sname
                        if 'LGZ' in s.sd[sname]['Created']:
                            if 'Renamed_from' in s.sd[sname]:
                                name=s.sd[sname]['Renamed_from']
                        s.log('Adding',name,'(component',cname,') to zoom list')
                        s.addzoom(name,'Rescue file')

    s.set_stage('Assoc check')
    # At this point Assoc should be correct except for sources where
    # e.g. a zoom file has over-ridden a de-blend. Fix those
    for source in tqdm(s.sd):
        if 'Deleted' in s.sd[source]:
            continue
        if len(s.sd[source]['Children'])>1:
            if 'Assoc' not in s.sd[source]:
                s.log('Fixing',source,'created by',s.sd[source]['Created'],'which has assoc unset but',len(s.sd[source]['Children']),'children')
            elif s.sd[source]['Assoc']!=len(s.sd[source]['Children']):
                s.log('Fixing',source,'created by',s.sd[source]['Created'],'which has assoc =',s.sd[source]['Assoc'],'but',len(s.sd[source]['Children']),'children')
                s.sd[source]['Assoc']=len(s.sd[source]['Children'])
        elif len(s.sd[source]['Children'])==1:
            if 'Assoc' in s.sd[source] and s.sd[source]['Assoc']!=0:
                s.log('Fixing',source,'created by',s.sd[source]['Created'],'which has assoc =',s.sd[source]['Assoc'],'but',len(s.sd[source]['Children']),'children')
                s.sd[source]['Assoc']=0
    return s

def generate_table(sd,columns,keep_deleted=False): 
    sources=sorted(sd.keys())
    headers=[]
    cols=[]
    for c,default in columns:
        print(c,'... ',end='')
        sys.stdout.flush()
        headers.append(c)
        column=[]
        for sn in sources:
            if 'Deleted' in sd[sn] and not keep_deleted:
                continue
                
            if c in sd[sn]:
                column.append(sd[sn][c])
            elif default is None:
                if 'Created' not in sd[sn]:
                    print(sd[sn])
                    raise RuntimeError('Source %s ** without creation info ** missing required key %s' % (sn,c))
                else:  
                    raise RuntimeError('Source %s (created by route %s) missing required key %s' % (sn,sd[sn]['Created'],c))
            else:
                column.append(default)
        cols.append(column)
    t=Table(cols,names=headers)
    print('Done!')
    return t

def sanity_check(s):
    errors=0
    for source in tqdm(s.sd,desc='Sources   '):
        if 'Deleted' not in s.sd[source]:
            if 'Children' not in s.sd[source]:
                s.log('Source',source,'has no children!')
                s.log(s.sd[source])
                errors+=1
            else:
                # add children/assoc integrity check
                nchildren=len(s.sd[source]['Children'])
                if 'Assoc' not in s.sd[source]:
                    if nchildren>1:
                        s.log('Source',source,'has unset assoc but more than one child')
                        errors+=1
                    else:
                        #fine
                        pass
                elif s.sd[source]['Assoc']==0 and nchildren!=1:
                    s.log('Source',source,'has assoc=0 but more than 1 child')
                    errors+=1
                elif s.sd[source]['Assoc']>0 and nchildren!=s.sd[source]['Assoc']:
                    s.log('Source',source,'has assoc', s.sd[source]['Assoc'],'but',nchildren,'children')
                    errors+=1
                              
                    
                for component in s.sd[source]['Children']:
                    if 'Deleted' in s.cd[component]:
                        s.log(source)
                        s.log('Deleted component ',component,'in child list')
                        s.log('Source created by',s.sd[source]['Created'],'with children',s.sd[source]['Children'])
                        if s.sd[source]['Created']=='Deblend':
                            s.log('(blend file was %s)' % s.sd[source]['Blend_file'])
                        s.log('Component created by',s.cd[component]['Created'],'and deleted for reason',s.cd[component]['Deleted'])
                        errors+=1
                    elif s.cd[component]['Parent']!=source:
                        s.log('Source',source,'has a child',component,'that reports a different parent',s.cd[component]['Parent'])
                        s.log('Source created by',s.sd[source]['Created'],'with children',s.sd[source]['Children']) 
                    else:
                        s.cd[component]['Checked']=True
                        for gaussian in s.cd[component]['Children']:
                            if 'Deleted' in s.gd[gaussian]:
                                s.log(source,component)
                                s.log('Deleted Gaussian',gaussian,'in child list')
                                errors+=1
                            else:
                                s.gd[gaussian]['Checked']=True

    for component in tqdm(s.cd,desc='Components'):
        if 'Deleted' not in s.cd[component] and 'Checked' not in s.cd[component]:
            if 'Parent' not in s.cd[component]:
                s.log('Component',component,'has no parent!!')
                s.log(s.cd[component])
            else:
                parent=s.cd[component]['Parent']
                s.log('Component',component,'neither deleted nor part of a source (alleged parent is',parent,' and creation route is',s.cd[component]['Created'],')')
                if parent in s.sd and 'Deleted' in s.sd[parent]:
                    s.log('Parent is deleted for reason',s.sd[parent]['Deleted'])
                elif parent not in s.sd:
                    s.log('Parent',parent,'does not exist')
                else:
                    s.log('Parent properties are',s.sd[parent])
            errors+=1

    for gaussian in tqdm(s.gd,desc='Gaussians '):
        if 'Deleted' not in s.gd[gaussian] and 'Checked' not in s.gd[gaussian]:
            s.log('Gaussian',gaussian,'neither deleted nor part of a source (alleged parent is',s.gd[gaussian]['Parent'],')')
            errors+=1

    print('Total errors',errors)

def write_table(outname,sd,columns,rename=None): 

    t=generate_table(sd,columns)
    if rename is not None:
        for old,new in rename:
            t[old].name=new
    t.write(outname,overwrite=True)

if __name__=='__main__':

    # find existing version
    g=sorted(glob.glob('sources-v*.fits'))
    if len(g)>0:
        infile=g[-1]
        version=infile.replace('sources-','').replace('.fits','')
    else:
        version='v0.1'

    # version specified on command line over-rides this
    if len(sys.argv)>1:
        for a in sys.argv[1:]:
            if a[0]=='v':
                version=a

    dir=os.getcwd()
    field=os.path.basename(dir)

    print('field is',field)
    s=make_structure(field,version=version)

    print('*** Sanity checking *** ')
    
    sanity_check(s)
    
    print('Constructing output table')
    columns=[('Source_Name',None),('RA',None),('DEC',None),('E_RA',np.nan),('E_DEC',np.nan),('Total_flux',None),('E_Total_flux',None),('Peak_flux',None),('E_Peak_flux',None),('S_Code',None),('Mosaic_ID',None),('Maj',np.nan),('Min',np.nan),('PA',np.nan),('E_Maj',np.nan),('E_Min',np.nan),('E_PA',np.nan),('DC_Maj',np.nan),('DC_Min',np.nan),('DC_PA',np.nan),('Isl_rms',np.nan),('FLAG_WORKFLOW',-1),('Prefilter',0),('NoID',0),('lr_fin',np.nan),('UID_L',""),('optRA',np.nan),('optDec',np.nan),('LGZ_Size',np.nan),('LGZ_Width',np.nan),('LGZ_PA',np.nan),('Assoc',0),('Assoc_Qual',np.nan),('Art_prob',np.nan),('Blend_prob',np.nan),('Hostbroken_prob',np.nan),('Imagemissing_prob',np.nan),('Zoom_prob',np.nan),('Created',None),('Position_from',None),('Renamed_from',"")]
    write_table('sources-'+version+'.fits',s.sd,columns)

    columns=[('Source_Name',None),('RA',None),('DEC',None),('E_RA',None),('E_DEC',None),('Total_flux',None),('E_Total_flux',None),('Peak_flux',None),('E_Peak_flux',None),('S_Code',None),('Mosaic_ID',None),('Maj',np.nan),('Min',np.nan),('PA',np.nan),('E_Maj',np.nan),('E_Min',np.nan),('E_PA',np.nan),('DC_Maj',np.nan),('DC_Min',np.nan),('DC_PA',np.nan),('Created',None),('Deblended_from',""),('Parent',None)]
    rename=[('Source_Name','Component_Name'),('Parent','Parent_Source')]
    write_table('components-'+version+'.fits',s.cd,columns,rename=rename)

    if len(sys.argv)>1 and 'save' in sys.argv:
        s.save('structure-'+version)
