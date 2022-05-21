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
from zoomprep_dr2 import tzi_sql
from astropy_healpix import HEALPix
from astropy import units as u

scale=3600.0 # units of catalogue are arcsec

class Interactive(object):
    def __init__(self,f,optra,optdec,ots,components,notcomponents):
        self.f=f
        self.optra=optra
        self.optdec=optdec
        self.orig_optra=optra
        self.orig_optdec=optdec
        self.ots=ots
        self.components=components
        self.notcomponents=notcomponents
        self.redraw()
        self.set_mode('m')
        self.size=np.nan
        self.blend=False
        
    def redraw(self):
        f.set_auto_refresh(False)
        # remove the layers if already drawn
        try:
            self.f.remove_layer('ID_marker')
        except:
            pass
        try:
            self.f.remove_layer('Component_ellipse')
        except:
            pass
        
        if not(np.isnan(self.optra)):
            self.f.show_markers(self.optra,self.optdec,marker='x',facecolor='magenta',edgecolor='magenta',linewidth=3,s=1500,zorder=300,layer='ID_marker')

        c=[]
        for r in self.ots:
            if r['Source_Name'] in self.components:
                c.append('green')
                #self.f.show_ellipses(r['RA'],r['DEC'],r['Maj']*2/scale,r['Min']*2/scale,angle=90+r['PA'],edgecolor='green',linewidth=3,zorder=101)
            elif r['Source_Name'] in self.notcomponents:
                c.append('cyan')
                #self.f.show_ellipses(r['RA'],r['DEC'],r['Maj']*2/scale,r['Min']*2/scale,angle=90+r['PA'],edgecolor='cyan',linewidth=3,zorder=101)
            else:
                #self.f.show_ellipses(r['RA'],r['DEC'],r['Maj']*2/scale,r['Min']*2/scale,angle=90+r['PA'],edgecolor='red',linewidth=3,zorder=101)
                c.append('red')
        self.f.show_ellipses(self.ots['RA'],self.ots['DEC'],self.ots['Maj']*2/scale,self.ots['Min']*2/scale,angle=90+self.ots['PA'],edgecolor=c,linewidth=3,zorder=101,layer='Component_ellipse') 
        f.refresh()
        
    def onclick(self,event):
        xp=event.xdata
        yp=event.ydata
        ra,dec=self.f.pixel2world(xp,yp)
        if event.button==2:
            # do something which depends on mode
            if self.mode=='o':
                self.optra=ra
                self.optdec=dec
                print('Optical ID at',ra,dec)
                self.redraw()
            elif self.mode=='m':
                sep=separation(ra,dec,self.ots['RA'],self.ots['DEC'])
                index=np.argmin(sep)
                name=self.ots[index]['Source_Name']
                if self.ots[index]['Source_Name'] in self.components:
                    self.components.remove(name)
                    print('removed component',name)
                else:
                    self.components.append(name)
                    print('added component',name)
                self.redraw()
            elif self.mode=='z':
                if not(np.isnan(self.oldra)):
                    self.size=separation(ra,dec,self.oldra,self.olddec)*3600
                    print('Size measured as',self.size,'arcsec')
                self.oldra=ra
                self.olddec=dec
            elif self.mode=='i':
                sep=separation(ra,dec,self.ots['RA'],self.ots['DEC'])
                index=np.argmin(sep)
                print(self.ots[index])
            else:
                raise NotImplementedError('Mode not recognised')

    def set_mode(self,mode):
        self.mode=mode
        if self.mode=='m':
            print('Selected mark/unmark source mode')
        elif self.mode=='o':
            print('Selected optical ID mode')
            self.optra=np.nan
            self.optdec=np.nan
            self.redraw()
        elif self.mode=='z':
            print('Selected size mode')
            self.oldra=np.nan
            self.olddec=np.nan
        elif self.mode=='i':
            print('Selected component inspector mode')
        else:
            raise NotImplementedError('Mode not recognised')
        
    def write(self,name):
        outfile=open(name+'.txt','w')
        if len(self.components)>0:
            outfile.write('## Components\n')
            for c in self.components:
                outfile.write(c+'\n')
            outfile.write('\n')
        if not(np.isnan(self.optra)):
            outfile.write('## OptID\n%f %f\n\n' % (self.optra,self.optdec))
        if not(np.isnan(self.size)):
            outfile.write('## Size\n%f\n\n' % (self.size))

        if self.blend:
            outfile.write('## Blend\n\n')
        outfile.close()

    def delete(self,name):
        outfile=open(name+'.txt','w')
        outfile.write('## Deleted\n')
        for c in self.components:
            outfile.write(c+'\n')
        outfile.close()

def load_opt(h):
    try:
        t=Table.read(cwd+'/optical_hpix_256/%i.fits' % h)
    except IOError:
        t=None
    return t
                
if __name__=='__main__':

    dir=os.getcwd()
    field=os.path.basename(dir)
    print('field is',field)
    table=field.replace('-','_')
    sql=tzi_sql(table)
    
    if not os.path.isdir('zoom'):
        os.mkdir('zoom') # store text files in here

    g=sorted(glob.glob('structure*-sources.pickle'))
    print('Using structure',g[-1])
    s=Source.load(g[-1].replace('-sources.pickle',''))

    print('Reading files')

    dynamic_opt=False
    if os.path.isdir('optical_hpix_256'):
        dynamic_opt=True
        print('Using dynamic loading of optical data')
        od={}
        hp = HEALPix(nside=256)
    else:
        optcat='optical.fits'
        gals=Table.read(optcat)
        if 'RA' in gals.colnames:
            gals['RA'].name='ra'
        if 'DEC' in gals.colnames:
            gals['DEC'].name='dec'

    columns=[('Source_Name',None),('RA',None),('DEC',None),('E_RA',None),('E_DEC',None),('Total_flux',None),('E_Total_flux',None),('Peak_flux',None),('E_Peak_flux',None),('S_Code',None),('Maj',np.nan),('Min',np.nan),('PA',np.nan),('E_Maj',np.nan),('E_Min',np.nan),('E_PA',np.nan),('DC_Maj',np.nan),('DC_Min',np.nan),('DC_PA',np.nan),('Created',None),('Parent',None)]
    ot=generate_table(s.cd,columns,keep_deleted=True)
    # large source table
    lt=ot[(ot['Total_flux']>0.005) & (ot['Maj']>10)]
    cwd=os.getcwd()
    os.chdir('zoom')

    mode='wise'
    quit=False
    ds9=False
    while not(quit):
        if ds9:
            os.system('killall ds9')
            ds9=False
        sourcename=sql.get_next()
        if sourcename is None:
            break
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
            print(sourcename,'already has a zoom file')
            record['complete']=1
            sql.set_object(sourcename,record)
            continue
        if len(s.get_comps(sourcename))==0:
            print(sourcename,'has no components!')
            record['complete']=1
            sql.set_object(sourcename,record)
            continue
        if os.system('grep -q '+sourcename+' *.txt')==0:
            print(sourcename,'already exists as a zoom file component!')
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
        while True:
            print('Scale factor is',scalefactor)
            if overlaygals:
                rsize=size*scalefactor/np.cos(dec*np.pi/180.0)
                dsize=size*scalefactor
                if dynamic_opt:
                    lhp=hp.lonlat_to_healpix(np.array([ra,ra-rsize,ra+rsize,ra-rsize,ra+rsize])*u.deg,np.array([dec,dec-dsize,dec-dsize,dec+dsize,dec+dsize])*u.deg)
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
                                tables.append(newopt)
                            else:
                                print('Healpix',p,'not found')
                    gals=vstack(tables)
                    if 'RA' in gals.colnames:
                        gals['RA'].name='ra'
                    if 'DEC' in gals.colnames:
                        gals['DEC'].name='dec'

                pwg=gals[(np.abs(gals['ra']-ra)<rsize) & (np.abs(gals['dec']-dec)<dsize)]
            else:
                pwg=None

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
                f=show_overlay(lhdu,whdu,ra,dec,size*scalefactor/initfactor,coords_color='red',coords_ra=s.sd[sourcename]['RA'],coords_dec=s.sd[sourcename]['DEC'],coords_lw=3,lw=2,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,block=False,interactive=False,drlimit=8000,peak=peak,plotpos=[(pwg,'+')])
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
                if name in s.cd and s.cd['Parent']!=sourcename:
                    notcomponents.append(name)
            
            # notcomponents=s.get_ncomps(sourcename) # list of all components in other sources -- this is too big
            print('components is',components)
            I=Interactive(f,ora,odec,ots,components,notcomponents)
            fig=plt.gcf()
            fig.canvas.mpl_connect('button_press_event', I.onclick)

            plt.ion()
            plt.show(block=False)
            stop=False
            while not(stop):
                print('Source %s, %i sources remaining' % (sourcename,sql.get_remaining()))
                print('(d)rop source, (m)ark components (default), mark an (o)ptical ID,\n   mark a si(z)e, set (b)lend, save (f)its, go to (n)ext, (Z)oom out, \n   (i)nspect, (T)oggle galaxy overlay, change opt (I)mage,\n   (F)avourite, (s)ave and continue or save and (q)uit?',end=' ')
                command=raw_input()
                replot=False
                if command=='s' or command=='q':
                    stop=True
                    I.write(sourcename)
                    if command=='q':
                        quit=True
                if command=='d':
                    stop=True
                    I.delete(sourcename)
                if command=='n':
                    stop=True
                if command in ['m','o','z','i']:
                    I.set_mode(command)
                if command=='b':
                    I.blend=True
                if command=='p':
                    continue
                if 'f' in command:
                    lhdu.writeto(sourcename+'.fits',overwrite=True)
                    print('Saved as',sourcename+'.fits')
                    subprocess.Popen(['/soft/bin/ds9',sourcename+'.fits'])
                    ds9=True
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
                    replot=True
                if 'I' in command:
                    if mode=='wise':
                        mode='optical'
                    else:
                        mode='wise'
                    replot=True
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
                s.set_components(sourcename,I.components)
                s.set_opt(sourcename,I.optra,I.optdec)
                s.set_size(sourcename,I.size)
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
