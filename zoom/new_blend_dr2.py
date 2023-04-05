#!/usr/bin/python

from __future__ import print_function
from astropy.table import Table,vstack
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import sys
import os
import glob
import subprocess
from subim import extract_subim
import matplotlib.pyplot as plt
from make_catalogue import Source,make_structure,generate_table
from separation import separation
from image_utils import find_bbox,get_mosaic_name
from get_fits import save_fits
from overlay import show_overlay
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors
from new_blend_prep_dr2 import blend_sql
from astropy_healpix import HEALPix
from astropy import units as u
import pyds9

scale=3600.0 # units of catalogue are arcsec

class Assoc(object):
    # This is the equivalent of GaussAssoc from blend_dr2_lgz.  It
    # carries around the members (and non-members) of the association
    # together with a flag to say whether they were originally
    # Gaussians or not, and an optical position.

    # A key point here is that we can't regenerate the table structure
    # that's used to make ots every time. That's not an issue for the
    # size finding code but it means that the Interactive redraw needs
    # to use the info stored in this object and not the ots
    # table. Thus we don't actually take ots as an argument here but
    # populate the initial list using 'components' and
    # 'notcomponents'.

    def __init__(self,s,components,notcomponents,optra,optdec):
        self.s=s
        self.optids={}
        self.assocs={}
        self.types={}
        self.hidden={}
        self.parents={}
        self.unchanged=True
        self.flagged=False
        for c in components:
            self.assocs[c]=1
            self.types[c]='C'
            self.hidden[c]=False
        for c in notcomponents:
            self.assocs[c]=-1 # not in any source
            self.types[c]='C'
            self.hidden[c]=False
        if optra is not None and not(np.isnan(optra)):
            self.optids[1]=(optra,optdec)
            
    def optpos(self,c):
        try:
            return self.optids[c]
        except KeyError:
            return None,None

    def comp_notcomp(self,c):
        comp=[]
        notcomp=[]
        for cc in self.assocs:
            if self.assocs[cc]==c:
                comp.append(cc)
            elif self.assocs[cc]!=0:
                notcomp.append(cc)
        return comp,notcomp
        
    def remove(self,name):
        self.unchanged=False
        if name in self.assocs:
            self.assocs[name]=0
        else:
            raise RunTimeError('Removing a component which is not in the list')

    def ismember(self,c,name):
        return (self.assocs[name]==c)

    def add(self,c,name):
        self.unchanged=False
        if name in self.assocs:
            self.assocs[name]=c
        else:
            raise RunTimeError('Adding a component which is not in the list')

    def set_optid(self,c,ra,dec):
        self.unchanged=False
        self.optids[c]=(ra,dec)

    def split(self,name):
        # Split a component into its component Gaussians
        if name not in self.assocs:
            raise RuntimeError('Asked to split '+name+' which is not a source')
        if self.types[name]=='G':
            print('Cannot split',name,'as it is a Gaussian already')
            return False
        children=self.s.cd[name]['Children']
        if len(children)<2:
            print('Cannot split',name,'as it has only',len(children),'Gaussian component(s)')
            return False
        # now add the Gaussian children to the list
        for c in children:
            self.assocs[c]=self.assocs[name]
            self.types[c]='G'
            self.hidden[c]=False
        self.hidden[name]=True # hide original component
        return True
        
    def write(self,name):
        outfile=open(name+'.txt','w')
        if self.flagged:
            outfile.write('## Flagged\n')
        elif self.unchanged:
            outfile.write('## Unchanged\n')
        else:
            outfile.write('## Components\n')
            for c in self.assocs:
                if self.assocs[c]!=-1 and not self.hidden[c]:
                    outfile.write("%s %s %i\n" % (c,self.types[c],self.assocs[c]))
            outfile.write('\n## Optical IDs\n')
            for k in self.optids:
                outfile.write("%i %f %f\n" % ((k,)+self.optids[k]))
        outfile.close()

    def delete(self,name):
        outfile=open(name+'.txt','w')
        outfile.write('## Deleted\n')
        for c in self.assocs:
            if self.assocs[c]!=-1:
                outfile.write(c+'\n')
        outfile.close()

    def dump(self):
        print('Components:')
        for i,cc in enumerate(self.assocs):
            print(i,cc,self.assocs[cc],self.types[cc],self.hidden[cc])
        print('Optical IDs')
        for k in self.optids:
            print(k,self.optids[k])

    def ots(self,c):
        # On the fly, make a table with the current components and their drawing attributes
        ra=[]
        dec=[]
        gmaj=[]
        gmin=[]
        pa=[]
        color=[]
        ls=[]
        sources=[]
        lw=[]
        for k in self.assocs:
            if self.hidden[k]: continue
            if self.types[k]=='C':
                ls.append('-')
                lw.append(3)
                db=self.s.cd
            elif self.types[k]=='G':
                ls.append('--')
                lw.append(5)
                db=self.s.gd
            sources.append(k)
            ra.append(db[k]['RA'])
            dec.append(db[k]['DEC'])
            gmaj.append(db[k]['Maj'])
            gmin.append(db[k]['Min'])
            pa.append(db[k]['PA'])
            if self.assocs[k]==c:
                color.append('green') # current source
            elif self.assocs[k]==0:
                color.append('red') # removed from the original source and unassigned
            elif self.assocs[k]==-1:
                color.append('white') # part of another source
            else:
                color.append('cyan') # a member of another source number
                
        return Table([sources,ra,dec,gmaj,gmin,pa,color,ls,lw],names=['Source_Name','RA','DEC','Maj','Min','PA','Color','LS','LW'])

class Interactive(object):
    def __init__(self,f,ga):
        self.f=f # AplPy figure
        self.ga=ga
        self.set_mode('m')
        self.size=np.nan
        self.blend=False
        self.c=1
        self.redraw()
        
    def redraw(self):
        self.f.set_auto_refresh(False)
        # remove the layers if already drawn
        try:
            self.f.remove_layer('ID_marker')
        except:
            pass
        try:
            self.f.remove_layer('Component_ellipse')
        except:
            pass

        ra=[]
        dec=[]
        ls=[]
        for k in self.ga.optids:
            if k==self.c:
                linestyle='-'
            else:
                linestyle='--'
            optra,optdec=self.ga.optpos(k)
            if ra is None:
                continue
            ra.append(optra)
            dec.append(optdec)
            ls.append(linestyle)
        if len(ra)!=0:
            self.f.show_markers(ra,dec,marker='x',facecolor='magenta',edgecolor='magenta',linewidth=3,linestyle=ls,s=1500,zorder=300,layer='ID_marker')

        ots=ga.ots(self.c)
        self.f.show_ellipses(ots['RA'],ots['DEC'],ots['Maj']*2/scale,ots['Min']*2/scale,angle=90+ots['PA'],edgecolor=ots['Color'],linewidth=ots['LW'],linestyle=ots['LS'],zorder=101,layer='Component_ellipse') 
        self.f.refresh()

    def plot_galaxies(self,gals):
        self.f.set_auto_refresh(False)
        try:
            self.f.remove_layer('Galaxies')
        except:
            pass
        if gals is not None:
            f.show_markers(gals['ra'],gals['dec'],marker='+',facecolor='white',edgecolor='white',linewidth=2,s=750,zorder=100,layer='Galaxies')
        self.f.refresh()
        
    def onclick(self,event):
        ots=self.ga.ots(self.c)
        xp=event.xdata
        yp=event.ydata
        ra,dec=self.f.pixel2world(xp,yp)
        if event.button==2:
            # do something which depends on mode
            if self.mode=='o':
                self.ga.set_optid(self.c,ra,dec)
                print('\nOptical ID at',ra,dec)
                self.mode='m'
                print('Switching back to mark components mode!')
                self.redraw()
            elif self.mode=='m':
                sep=separation(ra,dec,ots['RA'],ots['DEC'])
                index=np.argmin(sep)
                name=ots[index]['Source_Name']
                if self.ga.ismember(self.c,name):
                    self.ga.remove(name)
                    print('\nremoved component',name,'from source',self.c)
                else:
                    self.ga.add(self.c,name)
                    print('\nadded component',name,'to source',self.c)
                self.redraw()
            elif self.mode=='b':
                sep=separation(ra,dec,ots['RA'],ots['DEC'])
                index=np.argmin(sep)
                name=ots[index]['Source_Name']
                if self.ga.split(name):
                    self.redraw()
                self.mode='m'
                print('Switching back to mark components mode!')
            elif self.mode=='z':
                if not(np.isnan(self.oldra)):
                    self.size=separation(ra,dec,self.oldra,self.olddec)*3600
                    print('Size measured as',self.size,'arcsec')
                self.oldra=ra
                self.olddec=dec
            elif self.mode=='r':
                if np.isnan(self.c1ra):
                    self.c1ra=ra
                    self.c1dec=dec
                    print('Select second corner!')
                else:
                    # Swap so c1 is the lower RA,DEC
                    if self.c1ra>ra:
                        ra,self.c1ra=self.c1ra,ra
                    if self.c1dec>dec:
                        ra,self.c1dec=self.c1ra,dec
                    for r in ots:
                        name=r['Source_Name']
                        if r['RA']>self.c1ra and r['RA']<ra and r['DEC']>self.c1dec and r['DEC']<dec:
                            if not self.ga.ismember(self.c,name):
                                self.ga.add(self.c,name)
                                print('added component',name,'to source',self.c)
                    print('Switching back to mark components mode!')
                    self.redraw()
                    self.mode='m'
            elif self.mode=='i':
                sep=separation(ra,dec,ots['RA'],ots['DEC'])
                index=np.argmin(sep)
                print(ots[index])
            else:
                raise NotImplementedError('Mode not recognised')

    def set_mode(self,mode):
        self.mode=mode
        if self.mode=='m':
            print('Selected mark/unmark source mode')
        elif self.mode=='o':
            print('Selected optical ID mode')
        elif self.mode=='z':
            print('Selected size mode')
            self.oldra=np.nan
            self.olddec=np.nan
        elif self.mode=='i':
            print('Selected component inspector mode')
        elif self.mode=='b':
            print('Selected break source mode')
        elif self.mode=='r':
            print('Selected rectangle select mode. Choose first corner')
            self.c1ra=np.nan
            self.c2dec=np.nan
        else:
            raise NotImplementedError('Mode not recognised')
        
    '''
    #not used
    def write(self,name):
        outfile=open(name+'.txt','w')
        if len(self.assocs)>0:
            outfile.write('## Components\n')
            for c in self.assocs:
                outfile.write(c+'\n')
            outfile.write('\n')
        if not(np.isnan(self.optra)):
            outfile.write('## OptID\n%f %f\n\n' % (self.optra,self.optdec))
        if not(np.isnan(self.size)):
            outfile.write('## Size\n%f\n\n' % (self.size))

        if self.blend:
            outfile.write('## Blend\n\n')
        outfile.close()
    '''

def load_opt(h):
    try:
        t=Table.read(cwd+'/optical_hpix_256/%i.fits' % h)
    except IOError:
        t=None
    return t
                
if __name__=='__main__':

    dir=os.getcwd()
    field=os.path.basename(dir)
    print('Spawn a ds9')
    pyds9.ds9_xpans()
    ds9=pyds9.DS9('new_blend')
    print('field is',field)
    table=field.replace('-','_')
    sql=blend_sql(table)
    
    if not os.path.isdir('new_blend'):
        os.mkdir('new_blend') # store text files in here

    g=sorted(glob.glob('structure*-sources.pickle'))
    if len(g)==0:
        print('Looks like you might be in the wrong working directory?')
        raise RuntimeError('No source files found')
    print('Using structure',g[-1])
    s=Source.load(g[-1].replace('-sources.pickle',''))

    dynamic_opt=False
    if os.path.isdir('optical_hpix_256'):
        dynamic_opt=True
        print('Using dynamic loading of optical data')
        od={}
        hp = HEALPix(nside=256)
    else:
        print('Reading optical files')
        optcat='optical.fits'
        gals=Table.read(optcat)
        if 'RA' in gals.colnames:
            gals['RA'].name='ra'
        if 'DEC' in gals.colnames:
            gals['DEC'].name='dec'

    columns=[('Source_Name',None),('RA',None),('DEC',None),('E_RA',None),('E_DEC',None),('Total_flux',None),('E_Total_flux',None),('Peak_flux',None),('E_Peak_flux',None),('S_Code',None),('Maj',np.nan),('Min',np.nan),('PA',np.nan),('E_Maj',np.nan),('E_Min',np.nan),('E_PA',np.nan),('DC_Maj',np.nan),('DC_Min',np.nan),('DC_PA',np.nan),('Created',None),('Parent',None)]
    print('Generating the table')
    ot=generate_table(s.cd,columns,keep_deleted=True)
    # large source table
    lt=ot[(ot['Total_flux']>0.005) & (ot['Maj']>10)]
    cwd=os.getcwd()
    os.chdir('new_blend')
    mode='wise'
    quit=False
    while not(quit):
        sourcename=sql.get_next()
        if sourcename is None:
            break
        if sourcename not in s.sd:
            print('Source',sourcename,'does not exist!')
            record['complete']=1
            sql.set_object(sourcename,record)
            continue

        record=sql.get_object(sourcename)
        s.sd[sourcename]['lofarfile']=os.environ['IMAGEDIR']+'/'+record['lofarfile']
        if record['legacyfile'] is not None:
            s.sd[sourcename]['legacyfile']=os.environ['IMAGEDIR']+'/downloads/'+record['legacyfile']
        else:
            s.sd[sourcename]['legacyfile']=None 
        s.sd[sourcename]['wisefile']=os.environ['IMAGEDIR']+'/downloads/'+record['wisefile']
        if 'Deleted' in s.sd[sourcename]:
            print(sourcename,'has been deleted for reason',s.sd[sourcename]['Deleted'])
            record['complete']=1
            sql.set_object(sourcename,record)
            continue
        if os.path.isfile(sourcename+'.txt'):
            print(sourcename,'already has a new_blend file')
            record['complete']=1
            sql.set_object(sourcename,record)
            continue
        if len(s.get_comps(sourcename))==0:
            print(sourcename,'has no components!')
            record['complete']=1
            sql.set_object(sourcename,record)
            continue
            
        ra,dec=s.sd[sourcename]['RA'],s.sd[sourcename]['DEC']
        lofarfile=s.sd[sourcename]['lofarfile']
        wisefile=s.sd[sourcename]['wisefile']
        legacyfile=s.sd[sourcename]['legacyfile']
        
        marker_ra=None
        marker_dec=None

        title=sourcename

        components=s.get_comps(sourcename)
        print('Sourcename is',sourcename,'components is',components)
        mask=np.array([False]*len(ot))
        for c in components:
            mask|=(ot['Source_Name']==c)
        ctable=ot[mask]
        print(ctable)
        # resize the image to look for interesting neighbours
        iter=0
        dist=90
        while True:
            startra,startdec=ra,dec
            tcopy=lt
            tcopy['dist']=np.sqrt((np.cos(dec*np.pi/180.0)*(tcopy['RA']-ra))**2.0+(tcopy['DEC']-dec)**2.0)*3600.0
            tcopy=tcopy[tcopy['dist']<dist]
            if len(tcopy)==0:
                dist*=1.2
                continue

            # all original components need to be in there
            tcopy=vstack([ctable,tcopy])
            ra=np.mean(tcopy['RA'])
            dec=np.mean(tcopy['DEC'])

            if startra==ra and startdec==dec:
                break
            iter+=1
            if iter==10:
                break

        # now find the bounding box of the resulting collection
        ra,dec,size=find_bbox(tcopy)

        if np.isnan(size):
            ra,dec=s.sd[sourcename]['RA'],s.sd[sourcename]['DEC']
            size=60

        if size>200:
            size=200.0
            ra,dec=s.sd[sourcename]['RA'],s.sd[sourcename]['DEC']
        if size<60:
            size=60.0
        size=(int(0.5+size/10))*10

        size/=3600.0
        initfactor=4.0
        scalefactor=initfactor
        overlaygals=False
        ga=None
        while True:
            print('Scale factor is',scalefactor)

            ots=ot[separation(ra,dec,ot['RA'],ot['DEC'])<(size*scalefactor)]

            #pshdu=extract_subim(imagedir+'/downloads/'+psmaps[i],ra,dec,size*2,hduid=1)
            print('Lofarfile is',lofarfile)
            lhdu=extract_subim(lofarfile,ra,dec,size*scalefactor)
            if mode=='wise':
                whdu=extract_subim(wisefile,ra,dec,size*scalefactor)
            else:
                if legacyfile is None:
                    print('Legacy image does not exist, using WISE')
                    mode='wise'
                    whdu=extract_subim(wisefile,ra,dec,size*scalefactor)
                else:
                    whdu=extract_subim(legacyfile,ra,dec,size*scalefactor)
                    if np.all(np.isnan(whdu[0].data)):
                        print('Legacy image is blanked, using WISE')
                        mode='wise'
                        whdu=extract_subim(wisefile,ra,dec,size*scalefactor)

            try:
                peak==r['Peak_flux']/1000.0
            except:
                peak=None

            try:
                f=show_overlay(lhdu,whdu,ra,dec,size*scalefactor/initfactor,coords_color='red',coords_ra=s.sd[sourcename]['RA'],coords_dec=s.sd[sourcename]['DEC'],coords_lw=3,lw=2,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,block=False,interactive=False,drlimit=8000,peak=peak,rotate_north=False)
            except Exception as e:
                print('*** Error making image. Press RETURN to skip and continue ***')
                print(e)
                raw_input()
                stop=True
                break

            ora,odec=s.optid(sourcename)

            components=s.get_comps(sourcename)
            notcomponents=[]
            for r in ots:
                name=r['Source_Name']
                if name in s.cd and name not in components:
                    notcomponents.append(name)

            print('components is',components)
            # we'll have to think about how this works if we zoom out
            # and more components get into the ots list
            if ga is None:
                ga=Assoc(s,components,notcomponents,ora,odec)
            
            I=Interactive(f,ga)
            fig=plt.gcf()
            fig.canvas.mpl_connect('button_press_event', I.onclick)

            plt.ion()
            plt.show(block=False)
            stop=False
            while not(stop):
                print('Source %s, %i sources remaining' % (sourcename,sql.get_remaining()))
                print('Source number %i, display mode %s, operation mode %s' % (I.c,mode,I.mode))
                print('Number keys -- select new source number')
                print('(d)rop source, (m)ark components (default), mark (o)ptical ID,\n   (b)reak source, (r)ectangle select, save (f)its, (Z)oom out,\n   (i)nspect, (T)oggle galaxy overlay, change opt (I)mage,\n   (F)avourite, (s)ave and continue or save and (q)uit?',end=' ')
                command=raw_input()
                replot=False
                if command=='s' or command=='q':
                    stop=True
                    ga.write(sourcename)
                    if command=='q':
                        quit=True
                elif command>='0' and command<='9':
                    try:
                        I.c=int(command)
                    except ValueError:
                        print("Can't get number from input!")
                    #I.set_mode('O')
                    I.redraw()
                if command=='d':
                    stop=True
                    ga.delete(sourcename)
                if command=='n':
                    stop=True
                if command in ['m','o','z','i','b','r']:
                    I.set_mode(command)
                if command=='o' and I.c in ga.optids:
                    del(ga.optids[I.c])
                    ga.unchanged=False
                    I.redraw()
                if command=='p':
                    continue
                if 'f' in command:
                    ds9.set_pyfits(lhdu)
                    #lhdu.writeto(sourcename+'.fits',overwrite=True)
                    #print('Saved as',sourcename+'.fits')
                    #subprocess.Popen(['/soft/bin/ds9',sourcename+'.fits'])
                    #ds9=True
                if command=='F':
                    with open(os.environ['HOME']+'/favourites.txt','a') as outfile:
                        print('Adding',sourcename,'to favourites')
                        outfile.write(sourcename+'\n')
                if 'Z' in command:
                    zi=command.find('Z')
                    try:
                        scalefactor=int(command[zi+1:])
                    except:
                        scalefactor*=2
                    replot=True
                if 'T' in command:
                    overlaygals=not(overlaygals)
                    if overlaygals:
                        rsize=size*scalefactor/np.cos(dec*np.pi/180.0)
                        dsize=size*scalefactor
                        if dynamic_opt:
                            # do a cone search because looking at the corners may not work for large areas
                            lhp=hp.cone_search_lonlat(ra*u.deg,dec*u.deg,radius=size*u.deg)
                            lpix=sorted(list(set(lhp)))
                            print('Healpixes are',lpix)
                            tables=[]
                            for p in lpix:
                                if p in od:
                                    tables.append(od[p])
                                    print('Appending already loaded healpix',p)
                                else:
                                    newopt=load_opt(p)
                                    if newopt is not None:
                                        od[p]=newopt
                                        print('Dynamically loaded healpix',p)
                                        if 'RA' in newopt.colnames:
                                            newopt['RA'].name='ra'
                                        if 'DEC' in newopt.colnames:
                                            newopt['DEC'].name='dec'
                                        tables.append(newopt)
                                    else:
                                        print('Healpix',p,'not found')
                            gals=vstack(tables)

                        pwg=gals[(np.abs(gals['ra']-ra)<rsize) & (np.abs(gals['dec']-dec)<dsize)]
                    else:
                        pwg=None
                    I.plot_galaxies(pwg)
                if 'I' in command:
                    if mode=='wise':
                        mode='optical'
                    else:
                        mode='wise'
                    replot=True
                if command=='D':
                    ga.dump()
                    print(ga.ots(I.c))
                #else:
                #    print('Command not recognised!')
                if replot:
                    plt.close()
                    break

                    
            # out of inner loop
            if command=='d':
                s.delete_source(sourcename,'Marked as LGZ artefact')
            elif command!='n':
                # save what we've done here whether we are moving on or not; this preserves state if we change to a new zoom level half way through
                #s.set_components(sourcename,I.components)
                #s.set_opt(sourcename,I.optra,I.optdec)
                #s.set_size(sourcename,I.size)
                pass
            record['complete']=1
            sql.set_object(sourcename,record)
            
            if stop:
                break # out of outer while for scalefactor reset
            
        # options are:
        # drop source entirely
        # mark/unmark components
        # mark a size (for if components get this completely wrong)?
        # mark an optical ID (or just unmark existing if nothing selected)

        plt.close()
