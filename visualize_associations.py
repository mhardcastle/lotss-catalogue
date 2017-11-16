from astropy.table import Table
import shapely
from shapely.geometry import Polygon
from shapely.ops import cascaded_union
import numpy as np
import matplotlib.pyplot as plt
from descartes import PolygonPatch
import os

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

def separation(c_ra,c_dec,ra,dec):
    # all values in degrees
    return np.sqrt((np.cos(c_dec*np.pi/180.0)*(ra-c_ra))**2.0+(dec-c_dec)**2.0)

class visualize_compt(object):

    def plot(self,name,title=None):
        fig=plt.figure(1,figsize=(6,6))
        ax=fig.add_subplot(111)
        for p,c,a in zip(self.polys,self.colors,self.alphas):
            patch=PolygonPatch(p,color=c,alpha=a)
            ax.add_patch(patch)
        ax.set_aspect(1)
        plt.xlim(-self.field,self.field)
        plt.ylim(-self.field,self.field)
        if title: plt.title(title)
        plt.savefig(name+'.png')
        plt.clf()
        del(fig)

    def add_poly(self,poly,color,alpha):
        self.polys.append(poly)
        self.colors.append(color)
        self.alphas.append(alpha)
        
    def __init__(self,comps,ncomps):
        
        self.polys=[]
        self.colors=[]
        self.alphas=[]
        self.error=False

        # construct shape of old components
        ra=np.mean(comps['Comp_RA'])
        dec=np.mean(comps['Comp_DEC'])


        ellist=[]
        for r in comps:
            n_ra=r['Comp_RA']
            n_dec=r['Comp_DEC']
            x=3600*np.cos(dec*np.pi/180.0)*(ra-n_ra)
            y=3600*(n_dec-dec)
            newp=ellipse(x,y,r['Comp_Maj'],r['Comp_Min'],r['Comp_PA'])
            self.add_poly(newp,'red',0.3)
            ellist.append(newp)
        cp=cascaded_union(ellist)
            
        minx,miny,maxx,maxy=cp.bounds
        field=np.max([[maxx-minx],[maxy-miny]])
        self.field=field

        # find and remove small old components _not_ in the field --
        # catches compact sources in extended emission

        sep=3600.0*separation(ra,dec,ncomps['RA'],ncomps['DEC'])
        ncomps=ncomps[sep<field]
        for r in ncomps:
            n_ra=r['RA']
            n_dec=r['DEC']
            x=3600*np.cos(dec*np.pi/180.0)*(ra-n_ra)
            y=3600*(n_dec-dec)
            newp=ellipse(x,y,r['Maj'],r['Min'],r['PA'])
            self.add_poly(newp,'green',0.1)

t=Table.read('LOFAR_HBA_T1_DR1_catalog_v0.9.srl.fixed.fits')
oldt=Table.read('LOFAR_HBA_T1_DR1_catalog_v0.1.fits')
compt=Table.read('HETDEX-LGZ-comps-v0.3.fits')
sourcet=Table.read('HETDEX-LGZ-cat-v0.3-final.fits')

for i,r in enumerate(sourcet):
    if r['Assoc']==0:
        continue
    filter_com=compt['Source_Name']==r['Source_Name']
    comps=compt[filter_com]
    innew=[]
    for r2 in comps:
        name=r2['Comp_Name']
        innew.append(np.any(t['Source_Name']==name))
    if np.all(innew):
        print r['Source_Name'],'all matched in main cat'
        mask=np.array([True]*len(t))
        for r2 in comps:
            mask&=~(t['Source_Name']==r2['Comp_Name'])
        ncomps=t[mask]
    else:
        print r['Source_Name'],'partly or wholly unmatched with RA=',r['RA']
        mask=np.array([True]*len(oldt))
        for r2 in comps:
            mask&=~(oldt['Source_Name']==r2['Comp_Name'])
        ncomps=oldt[mask]
        assert(len(ncomps)+len(comps)==len(oldt))

    o=visualize_compt(comps,ncomps)
    name=r['Source_Name']
    o.plot(name)
    regionfile=name+'.png'
    if np.all(innew):
        dir='/data/lofar/mjh/hetdex_v4/extended_new/'
    else:
        dir='/data/lofar/mjh/hetdex_v3/lgz/subset_with_id_dir/'
    for r2 in comps:
        compfile=dir+r2['Comp_Name']+'_PS.png'
        print compfile
        if os.path.isfile(compfile):
            break
    else:
        print 'Argh argh argh'
        continue

    outfile=name+'_join.png'
    print regionfile,compfile
    os.system('convert +append '+regionfile+' '+compfile+' '+outfile)
