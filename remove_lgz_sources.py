from astropy.table import Table
import shapely
from shapely.geometry import Polygon
from shapely.ops import cascaded_union
import numpy as np
import matplotlib.pyplot as plt
from descartes import PolygonPatch
from copy import deepcopy

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

class find_overlap(object):

    def plot(self,name,title=None):
        fig=plt.figure(1,figsize=(6,6))
        ax=fig.add_subplot(111)
        for p,c,a in zip(self.polys,self.colors,self.alphas):
            try:
                patch=PolygonPatch(p,color=c,alpha=a)
                ax.add_patch(patch)
            except ValueError:
                print 'Some problem with polygon'
                print p
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
        
    def __init__(self,comps,ncomps,t):
        
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
            if r['Maj']>24:
                continue
            n_ra=r['RA']
            n_dec=r['DEC']
            x=3600*np.cos(dec*np.pi/180.0)*(ra-n_ra)
            y=3600*(n_dec-dec)
            newp=ellipse(x,y,r['Maj'],r['Min'],r['PA'])
            self.add_poly(newp,'green',0.1)
            cp=cp.difference(newp)

        self.add_poly(cp,'red',0.5)
        
        # compute nearby new components
        sep=3600.0*separation(ra,dec,t['RA'],t['DEC'])
        t=t[sep<field]
        print '        .... total nearby components',len(t)
        # set up ellipses for association. units are arcsec

        result=[]
        for r in t:
            n_ra=r['RA']
            n_dec=r['DEC']
            x=3600*np.cos(dec*np.pi/180.0)*(ra-n_ra)
            y=3600*(n_dec-dec)
            newp=ellipse(x,y,r['Maj'],r['Min'],r['PA'])
            error=False
            try:
                inter=cp.intersection(newp)
            except shapely.geos.TopologicalError:
                print 'Error!',x,y,r['Maj'],r['Min'],r['PA']
                error=True    
            self.add_poly(newp,'blue',0.2)
            if error:
                result.append(True)
                self.error=True
            elif inter.area>(0.4*newp.area):
                # do new sources overlap with old
                result.append(True)
                self.add_poly(newp,'yellow',0.2)
            else:
                # is old source subsumed by one big new source
                inter=newp.intersection(cp)
                if inter.area>(0.9*cp.area):
                    result.append(True)
                    self.add_poly(newp,'yellow',0.2)
                else:
                    result.append(False)
        self.filtered=t[result]

t=Table.read('LOFAR_HBA_T1_DR1_catalog_v0.9.srl.fixed.fits')
oldt=Table.read('LOFAR_HBA_T1_DR1_catalog_v0.1.fits')
compt=Table.read('HETDEX-LGZ-comps-v0.3.fits')
sourcet=Table.read('HETDEX-LGZ-cat-v0.3-final.fits')
sourcet['Component_flux']=np.nan
matched=0
unmatched=0
bad=0
remove=[]
drops=[]
for i,r in enumerate(sourcet):
    drop=False
    old_com=compt['Source_Name']==r['Source_Name']
    comps=compt[old_com]
    innew=[]
    for r2 in comps:
        name=r2['Comp_Name']
        innew.append(np.any(t['Source_Name']==name))
    if np.all(innew):
        print r['Source_Name'],'all matched in main cat'
        matched+=1
        for r2 in comps:
            remove.append(r2['Comp_Name'])
    else:
        print r['Source_Name'],'partly or wholly unmatched with RA=',r['RA']
        # now search for components overlapping with this one in the old table
        mask=np.array([True]*len(oldt))
        for r2 in comps:
            mask&=~(oldt['Source_Name']==r2['Comp_Name'])
        ncomps=oldt[mask]
        assert(len(ncomps)+len(comps)==len(oldt))
        o=find_overlap(comps,ncomps,t)
        nt=o.filtered
        print '        .... overlapping components are',len(nt)
        flux=0
        for r2 in nt:
            flux+=r2['Total_flux']
        print '        .... sanity check, source flux =',r['Flux'],'comps total flux =',flux
        sourcet[i]['Component_flux']=flux
        if flux==0:
            bad+=1
            drop=True
            #o.plot()
        elif flux/r['Flux']<0.5 or flux/r['Flux']>2.0:
            bad+=1
            if flux/r['Flux']>3.0:
                drop=True
        unmatched+=1
        if drop: o.plot(r['Source_Name'],'Catalogue flux %f New component flux %f' % (r['Flux'],flux))
        else: remove.append(r2['Source_Name'])
    drops.append(drop)

print matched,unmatched,bad

stf=sourcet[drops]
stf.write('drop.fits',overwrite=True)

removenames=open('remove.txt','w')
for r in remove:
    removenames.write(r+'\n')
removenames.close()
