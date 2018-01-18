# Code to process LGZ output with the too-zoomed-in (or other)
# datasets. Reads in the output from Jude's code plus the per-source
# data that are generated by zoom.py. Then, if sources in Jude's list
# are unaltered, we pass them through unchanged, but if they are
# affected by the per-source files, we add or remove lines from the
# table as appropriate.

from astropy.table import Table
import numpy as np
import os
from copy import deepcopy
from astropy.coordinates import SkyCoord
from separation import separation
import astropy.units as u
import glob

import shapely
from shapely.geometry import Polygon
from shapely.ops import cascaded_union

from source_handler import Source,parsefile

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

if __name__=='__main__':

    scale=3600.0 # scaling factor for region sizes

    do_optical=True # Set to false for testing purposes only

    print 'Reading tables, please wait...'
    pyb=Table.read('../LOFAR_HBA_T1_DR1_catalog_v0.95_masked.srl.fits')
    lgz=Table.read('HETDEX-LGZ-cat-v0.6-filtered.fits')
    lgz['New_size']=np.nan
    lgz['New_width']=np.nan
    lgz['New_PA']=np.nan
    if do_optical:
        ocat=Table.read('/data/lofar/mjh/hetdex_ps1_allwise_photoz_v0.2.fits')
        opcat=SkyCoord(ocat['ra']*u.deg, ocat['dec']*u.deg)

    ss=Source()

    for i,r in enumerate(lgz):
        k=r['Source_Name']
        ss.set_components(k,[])
        ss.set_lgz_number(k,i)

    # components.txt contains components for straight lgz outputs

    clines=open('components.txt').readlines()
    for c in clines:
        bits=c.split()
        ss.add(bits[0],bits[1].rstrip())
        ss.mdict[bits[0]]=1

    ss.reset_changes() # ignore all 'changes' from initialization

    # zooms table includes overrides for all sources, including ones that may not have come from lgz

    dir='/data/lofar/mjh/hetdex_v4/zoom/'
    g=glob.glob(dir+'ILT*.txt')
    for f in g:
        if 'table-list' in f:
            continue
        n=f.split('/')[-1].replace('.txt','')
        print 'Parsing file for source',n
        parsefile(n,ss,dir=dir)
        ss.mdict[n]=2

    olgz=lgz[:0].copy()

    remove=open('lgz_components.txt','w')

    for k in ss.cdict:
        if len(ss.cdict[k])==0:
            print 'Dropping source',k,'as it has no components'
            r=None
        else:
            # Get source components, whether it's changed or not
            pfilter=np.array([False]*len(pyb))
            cnames=ss.get_comps(k)
            for n in cnames:
                pfilter|=(pyb['Source_Name']==n)
            clist=pyb[pfilter]
            assert(len(cnames)==len(clist))
            ms=Make_Shape(clist)

            if not(ss.changed_dict[k]):
                print 'Source',k,'has not changed'
                # source has not changed, which means it can be copied from the corresponding entry in the lgz table
                r=lgz[ss.idict[k]]
            else:
                # source has changed, or was only ever present in zooms text
                print 'Source',k,'is flagged as having changed'
                r=lgz[0]
                # Clear some fields to default values
                for j in ['Art_prob','Blend_prob','Zoom_prob','Hostbroken_prob','ID_Qual','Assoc_Qual','Assoc', 'Compoverlap']:
                    r[j]=0
                r['OptID_Name']='None'
                r['optRA']=np.nan
                r['optDec']=np.nan
                r['Assoc_Qual']=1 # these should come from zoom files and are therefore perfect
                if k in ss.blends:
                    r['Blend_prob']=1

            r['New_size']=ms.length()
            r['New_width']=ms.width()
            r['New_PA']=ms.pa()
            # WHETHER THE SOURCE HAS CHANGED OR NOT,
            # recreate all of the source's radio info from its component list
            # this ensures consistency
            try:
                ora,odec=ss.odict[k]
            except KeyError:
                ora=None
                odec=None
            if ora is not None:
                print '      Updated optical position',ora,odec
                # later we check for source near this position
            try:
                size=ss.sdict[k]
            except KeyError:
                size=None

            if size is not None:
                print '      Updated size',size

            r['Assoc']=len(clist)
            if len(clist)==1:
                # only one component, so we can use its properties
                c=clist[0]
                for key in ['RA','DEC','Total_flux','E_Total_flux','Peak_flux','E_Peak_flux','E_RA','E_DEC','Isl_rms','S_Code','Mosaic_ID','Source_Name']:
                    r[key]=c[key]
                r['Size']=c['DC_Maj']
                if size is not None:
                    r['Size']=size
                    if size>r['New_size']:
                        r['New_size']=size
            else:
                tfluxsum=np.sum(clist['Total_flux'])
                ra=np.sum(clist['RA']*clist['Total_flux'])/tfluxsum
                dec=np.sum(clist['DEC']*clist['Total_flux'])/tfluxsum
                sname=sourcename(ra,dec)
                print '      New sourcename is',sname
                r['RA']=ra
                r['DEC']=dec
                r['Source_Name']=sname
                r['E_RA']=np.sqrt(np.mean(clist['E_RA']**2.0))
                r['E_DEC']=np.sqrt(np.mean(clist['E_DEC']**2.0))
                r['Source_Name']=sname
                r['Total_flux']=np.sum(clist['Total_flux'])
                r['E_Total_flux']=np.sqrt(np.sum(clist['E_Total_flux']**2.0))
                maxpk=np.argmax(clist['Peak_flux'])
                r['Peak_flux']=clist[maxpk]['Peak_flux']
                r['E_Peak_flux']=clist[maxpk]['E_Peak_flux']
                r['S_Code']='M'
                r['Isl_rms']=np.mean(clist['Isl_rms'])
                r['Mosaic_ID']=clist[maxpk]['Mosaic_ID']
                seps=[]
                for c in clist:
                    seps.append(separation(c['RA'],c['DEC'],clist['RA'],clist['DEC']))
                maxsep=np.max(seps)*scale
                maxsize=np.max(clist['Maj'])
                maxsize=max((maxsep,maxsize))
                if size is not None:
                    if size>maxsize:
                        maxsize=size
                    if size>r['New_size']:
                        r['New_size']=size

                print 'sizes:',maxsep,maxsize
                r['Size']=maxsize

            r['Zoom_prob']=0 # since these should all have been resolved
            if ora is not None and do_optical:
                # check opt position
                sep=separation(ora,odec,r['optRA'],r['optDec'])
                if np.isnan(sep) or sep>1.0:
                    print '         Separation is',sep,'updating with new opt pos'
                    c=SkyCoord(ora*u.deg,odec*u.deg)
                    idx, d2d, d3d = c.match_to_catalog_sky(opcat)
                    sep=d2d.value*scale
                    if sep>12:
                        print '     Bad optical position! sep = ',sep
                    else:
                        idx=int(idx)
                        ora=ocat[idx]['ra']
                        odec=ocat[idx]['dec']
                        if ocat[idx]['AllWISE']!='N/A':
                            name='AllWISE'+ocat[idx]['AllWISE']
                        else:
                            name='PSO %s' % ocat[idx]['objID']
                        qual=0.667+0.333*(12.0-sep)/12.0

                        r['optRA']=ora
                        r['optDec']=odec
                        r['OptID_Name']=name
                        r['ID_Qual']=qual

        if r is not None:
            olgz.add_row(r)

            comps=ss.get_comps(k)
            for j in comps:
                remove.write('%s %s %i\n' % (j,r['Source_Name'],ss.mdict[k]))

    olgz.write('HETDEX-LGZ-cat-v0.9-filtered-zooms.fits',overwrite=True)

    remove.close()
