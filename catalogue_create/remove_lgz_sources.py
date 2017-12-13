from astropy.table import Table
import shapely
from shapely.geometry import Polygon
from shapely.ops import cascaded_union
import numpy as np
import matplotlib.pyplot as plt
from descartes import PolygonPatch
from copy import deepcopy
from separation import separation

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
compt=Table.read('HETDEX-LGZ-comps-lgzv1.fits')
sourcet=Table.read('HETDEX-LGZ-cat-lgzv1.fits')
for name in ['Component_flux','E_RA','E_DEC','Peak_flux','E_Peak_flux','E_Total_flux','Isl_rms']:
    sourcet[name]=np.nan
sourcet['S_Code']=""
sourcet['Mosaic_ID']="          "
sourcet['Dec'].name='DEC'
sourcet['Flux'].name='Total_flux'

matched=0
unmatched=0
matched_in=0
matched_out=0
unmatched_in=0
unmatched_out=0
bad=0
remove=[]
rsources=[]
drops=[]
for i,r in enumerate(sourcet):
    drop=False
    select_com=compt['Source_Name']==r['Source_Name']
    comps=compt[select_com]
    innew=[]
    filter=np.array([False]*len(t))
    for r2 in comps:
        name=r2['Comp_Name']
        filter|=(t['Source_Name']==name)
    if np.sum(filter)==len(comps):
        print r['Source_Name'],'all matched in main cat'
        matched+=1
        matched_in+=len(comps)
        complist=t[filter]
        for r2 in comps:
            rsources.append(r['Source_Name'])
            remove.append(r2['Comp_Name'])
            matched_out+=1
    else:
        print r['Source_Name'],'partly or wholly unmatched with RA=',r['RA']
        # now search for components overlapping with this one in the old table
        unmatched_in+=len(comps)
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
        sourceflux=r['Total_flux']
        print '        .... sanity check, source flux =',sourceflux,'comps total flux =',flux
        sourcet[i]['Component_flux']=flux
        if flux==0:
            bad+=1
            drop=True
            #o.plot()
        elif flux/sourceflux<0.5 or flux/sourceflux>2.0:
            bad+=1
            if flux/sourceflux>3.0:
                drop=True
        unmatched+=1
        complist=nt
        if drop:
            o.plot(r['Source_Name'],'Catalogue flux %f New component flux %f' % (sourceflux,flux))
        if drop:
            if r['Compoverlap']>0 or r['Art_prob']>=0.5 or r['Blend_prob']>=0.5 or r['Hostbroken_prob']>0.5 or r['Zoom_prob']>=0.5:
                # these should be removed from the flowchart anyway because we need to fix them up later in various ways, so let them through now and sort them out later
                drop=False
        if not(drop):
            for r2 in nt:
                unmatched_out+=1
                rsources.append(r['Source_Name'])
                remove.append(r2['Source_Name'])
    drops.append(drop)
    # fill in missing columns
    # Assume total_flux is correct
    # complist has whatever we know about the pybdsf components in it
    if len(complist)>0:
        # errors on RA, DEC are errors on the mean
        sourcet[i]['E_RA']=np.sqrt(np.sum(complist['E_RA']**2.0))/len(complist)
        sourcet[i]['E_DEC']=np.sqrt(np.sum(complist['E_DEC']**2.0))/len(complist)
        # total flux error is error on the sum
        sourcet[i]['E_Total_flux']=np.sqrt(np.sum(complist['E_Total_flux']**2.0))
        # peak flux and error from brightest component
        maxpk=np.argmax(complist['Peak_flux'])
        sourcet[i]['Peak_flux']=complist[maxpk]['Peak_flux']
        sourcet[i]['E_Peak_flux']=complist[maxpk]['E_Peak_flux']

        # S_Code is determined from how many components we had
        if len(complist)==1:
            sourcet[i]['S_Code']=complist[0]['S_Code']
        else:
            sourcet[i]['S_Code']='M'

        # isl_rms estimated from all the islands
        sourcet[i]['Isl_rms']=np.mean(complist['Isl_rms'])

        # mosaic id from brightest component again
        sourcet[i]['Mosaic_ID']=complist[maxpk]['Mosaic_ID']
        
print 'Total matched',matched,'Total unmatched',unmatched,'Bad fluxes',bad
print 'Matched in',matched_in,'Matched out',matched_out
print 'Unmatched in',unmatched_in,'Unmatched out',unmatched_out

drops=np.array(drops)
stf=sourcet[drops]
stf.write('drop.fits',overwrite=True)

sourcet.remove_column('Component_flux')
stf=sourcet[~drops]
stf.write('HETDEX-LGZ-cat-v0.6-filtered.fits',overwrite=True)

removenames=open('remove.txt','w')
for r in remove:
    removenames.write(r+'\n')
removenames.close()

components=open('components.txt','w')
for s,r in zip(rsources,remove):
    components.write(s+' '+r+'\n')
components.close()
