from __future__ import print_function
from astropy.table import Table,vstack
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import sys
import os
import glob
from subim import extract_subim
import matplotlib.pyplot as plt
from separation import separation
from image_utils import find_bbox,get_mosaic_name
from overlay import show_overlay
from download_image_files import LofarMaps,get_legacy,get_first,get_wise
from find_wise import WISE
from builtins import input

scale=3600.0 # units of catalogue are arcsec

class GaussAssoc(object):
    # This class carries around the division of the source into components
    def __init__(self,ots,optra,optdec):
        self.ots=ots
        self.optids={}
        self.assocs={}
        self.unchanged=True
        self.flagged=False
        self.components=[r['Gaus_id'] for r in ots]
        for c in self.components:
            self.assocs[c]=1
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
        for cc in self.components:
            if self.assocs[cc]==c:
                comp.append(cc)
            elif self.assocs[cc]!=0:
                notcomp.append(cc)
        return comp,notcomp
        
    def remove(self,name):
        self.unchanged=False
        if name in self.components:
            self.assocs[name]=0
        else:
            raise RunTimeError('Removing a component which is not in the list')

    def ismember(self,c,name):
        return (self.assocs[name]==c)

    def add(self,c,name):
        self.unchanged=False
        if name in self.components:
            self.assocs[name]=c
        else:
            raise RunTimeError('Adding a component which is not in the list')

    def set_optid(self,c,ra,dec):
        self.unchanged=False
        self.optids[c]=(ra,dec)
        
    def write(self,name):
        outfile=open('blend/'+name+'.txt','w')
        if self.flagged:
            outfile.write('## Flagged\n')
        elif self.unchanged:
            outfile.write('## Unchanged\n')
        else:
            outfile.write('## Components\n')
            for c in self.components:
                outfile.write("%i %i\n" % (c,self.assocs[c]))
            outfile.write('\n## Optical IDs\n')
            for k in self.optids:
                outfile.write("%i %f %f\n" % ((k,)+self.optids[k]))
        outfile.close()

    def dump(self):
        print('Components:')
        for i,cc in enumerate(self.components):
            print(i,cc,self.assocs[cc],self.ots[i]['Total_flux'])
        print('Optical IDs')
        for k in self.optids:
            print(k,self.optids[k])
        
class Interactive(object):
    def __init__(self,f,ga):
        self.f=f
        self.ga=ga
        self.ots=ga.ots
        self.set_mode('m')
        self.c=1
        self.redraw()
        
    def redraw(self):
        f.set_auto_refresh(False)
        optra,optdec=self.ga.optpos(self.c)
        if optra is not None:
            self.f.show_markers(optra,optdec,marker='x',facecolor='magenta',edgecolor='magenta',linewidth=3,s=1500,zorder=300,layer='ID_marker')
        else:
            try:
                self.f.remove_layer('ID_marker')
            except:
                pass
        components,notcomponents=self.ga.comp_notcomp(self.c)
        c=[]
        for r in self.ots:
            if r['Gaus_id'] in components and self.c!=0:
                c.append('green')
                #self.f.show_ellipses(r['RA'],r['DEC'],r['Maj']*2/scale,r['Min']*2/scale,angle=90+r['PA'],edgecolor='green',linewidth=3,zorder=101)
            elif r['Gaus_id'] in notcomponents:
                c.append('cyan')
                #self.f.show_ellipses(r['RA'],r['DEC'],r['Maj']*2/scale,r['Min']*2/scale,angle=90+r['PA'],edgecolor='cyan',linewidth=3,zorder=101)
            else:
                #self.f.show_ellipses(r['RA'],r['DEC'],r['Maj']*2/scale,r['Min']*2/scale,angle=90+r['PA'],edgecolor='red',linewidth=3,zorder=101)
                c.append('red')
        self.f.show_ellipses(self.ots['RA'],self.ots['DEC'],self.ots['Maj']*2/scale,self.ots['Min']*2/scale,angle=90+self.ots['PA'],edgecolor=c,linewidth=3,zorder=101,layer='Component_ellipse') 
        f.refresh()

    def clickon(self,ra,dec):
        # do something which depends on mode
        if self.mode=='o' or self.mode=='O':
            self.ga.set_optid(self.c,ra,dec)
            print('\nOptical ID at',ra,dec)
            if self.mode=='o': self.redraw()
        if self.mode=='m' or self.mode=='O':
            sep=separation(ra,dec,self.ots['RA'],self.ots['DEC'])
            index=np.argmin(sep)
            name=self.ots[index]['Gaus_id']
            if self.ga.ismember(self.c,name):
                self.ga.remove(name)
                print('\nremoved component',name,'from source',self.c)
            else:
                self.ga.add(self.c,name)
                print('\nadded component',name,'to source',self.c)
            self.redraw()
            if self.mode=='O':
                self.set_mode('m')

        
    def onclick(self,event):
        xp=event.xdata
        yp=event.ydata
        ra,dec=self.f.pixel2world(xp,yp)
        if event.button==2:
            self.clickon(ra,dec)

    def set_mode(self,mode):
        self.mode=mode
        if self.mode=='m':
            print('Selected mark/unmark source mode')
        elif self.mode=='o':
            print('Selected optical ID mode')
        elif self.mode=='O':
            print('Selected select-and-ID mode')
        else:
            raise NotImplementedError('Mode not recognised')

if __name__=='__main__':

    dir=os.getcwd()
    field=os.path.basename(dir)
    print('field is',field)
    if not os.path.isdir('blend'):
        os.mkdir('blend') # store text files in here

    t=Table.read('LGZ-cat.fits')
    filt=t['Blend_prob']>0.5
    filt|=(t['Badclick']>=2) & (t['Zoom_prob']<=0.5)
    filt|=t['OptID_Name']=='Mult'
    t=t[filt]
    print('About to deblend',len(t),'sources:',len(glob.glob('blend/*.txt')),'blend files already exist')
    mt=Table.read('LGZ-multiple.fits')
    print('Reading data...')
    sourcecat='source_lr.fits'
    optcat='optical.fits'
    opt=Table.read(optcat)
    idname='ID'
    st=Table.read(sourcecat)
    gt=Table.read('edited_gaussians.fits') # has parent name and source name
    gt['ra']=np.where(gt['ra']>360,np.nan,gt['ra'])
    gt['dec']=np.where(gt['dec']>90,np.nan,gt['dec'])
    # now download any images we need!

    lofarfiles={}
    legacyfiles={}
    wisefiles={}
    w=WISE()
    if not os.path.isfile('blend-list.txt'):
        outfile=open('blend-list.txt','w')
        # now work in downloads dir
        wd=os.getcwd() # the zoom directory
        lm=LofarMaps(stay_in_imagedir=True)
        os.chdir('downloads')

        print('Making images')
        for r in t:
            sourcename=r['Source_Name']
            ra=r['RA']
            dec=r['Dec']
            lofarfiles[sourcename]=os.environ['IMAGEDIR']+'/'+lm.find(ra,dec)
            legacyfiles[sourcename]=os.environ['IMAGEDIR']+'/downloads/'+get_legacy(ra,dec,bands='zrg')
            wisename=w.find_pos(ra,dec)
            if wisename is None:
                wisename=get_wise(ra,dec,1)
            wisefiles[sourcename]=os.environ['IMAGEDIR']+'/downloads/'+wisename
            print(sourcename,lofarfiles[sourcename],legacyfiles[sourcename],wisefiles[sourcename],file=outfile)
        outfile.close()
        # back to zoom dir
        os.chdir(wd)
    else:
        lines=open('blend-list.txt').readlines()
        for l in lines:
            bits=l.rstrip().split()
            sourcename=bits[0]
            lofarfiles[sourcename]=bits[1]
            legacyfiles[sourcename]=bits[2]
            wisefiles[sourcename]=bits[3]
        
    gals=opt
    if 'RA' in gals.colnames:
        gals['RA'].name='ra'
    if 'DEC' in gals.colnames:
        gals['DEC'].name='dec'
        
    # read lists
    display_mode='wise'

    for i,r in enumerate(t):
        name=r['Source_Name']
        lofarfile=lofarfiles[name]
        wisefile=wisefiles[name]
        legacyfile=legacyfiles[name]
        print('Optid name is',r['OptID_Name'])
        # work out whether there are sources selected in LGZ
        if r['OptID_Name']=="None":
            selected=None
        elif r['OptID_Name']=="Mult":
            mts=mt[mt['Source_Name']==name]
            selected=[row['Opt_id'] for row in mts]
        else:
            selected=[r['OptID_Name']]
        if selected is not None:
            gzlist=[]
            for s in selected:
                gzlist.append(np.argmax(gals[idname]==s))
            print(gzlist)
            gzgals=gals[gzlist]
        else:
            gzgals=None
        print('Selected list is',selected)
            
        sts=st[st['Source_Name']==name]
        if len(sts)==0:
            print(name,'does not exist in source table, skipping')
            continue
        r2=sts[0] # corresponding row in source table
        if os.path.isfile('blend/'+name+'.txt'):
            print(name,'already has a blend file')
            continue
        print(r2)
        rgt=gt[gt['Parent_Name']==name]
        print('Contains',len(rgt),'Gaussians')
        mlt=rgt[~np.isnan(rgt['ra'])]
        mlt['ra']=mlt['ra']
        mlt['dec']=mlt['dec']
        title=name+' (%s)' % display_mode
        ra,dec,size=find_bbox(rgt,scale=scale)
        size*=1.5
        size/=3600.0
        print('Lofarfile is',lofarfile)
        lhdu=extract_subim(lofarfile,ra,dec,size*2)
        try:
            peak==r2['Peak_flux']/1000.0
        except:
            peak=None
        pwg=gals[(np.abs(gals['ra']-ra)<size) & (np.abs(gals['dec']-dec)<size)]

        if gzgals is not None:
            # pick a LGZ ID over the LR one
            ora=gzgals[0]['ra']
            odec=gzgals[0]['dec']
        else:
            ora=r2['ra']
            odec=r2['dec']
        ga=GaussAssoc(rgt,ora,odec)

        ngt=gt[(np.abs(gt['RA']-ra)<size) & (np.abs(gt['DEC']-dec)<size)]
        ngt=ngt[ngt['Parent_Name']!=name]
        
        stop=False
        while not(stop):
            if display_mode=='wise':
                ohdu=extract_subim(wisefile,ra,dec,size)
            else:
                ohdu=extract_subim(legacyfile,ra,dec,size)
            f=show_overlay(lhdu,ohdu,ra,dec,size,coords_color='red',coords_ra=r['RA'],coords_dec=r['Dec'],coords_lw=3,lw=2,no_labels=True,title=title,block=False,interactive=False,drlimit=8000,peak=peak,plotpos=[(pwg,'+'),(mlt,'h','red','none'),(gzgals,'D','magenta','none')])
            if len(ngt)>0:
                f.show_ellipses(ngt['RA'],ngt['DEC'],ngt['Maj']*2/scale,ngt['Min']*2/scale,angle=90+ngt['PA'],edgecolor='white',linewidth=2,zorder=10,layer='BG_ellipse')

            I=Interactive(f,ga)
            fig=plt.gcf()
            fig.canvas.mpl_connect('button_press_event', I.onclick)

            plt.ion()
            plt.show(block=False)

            while not(stop):
                print('Object number %i of %i' % (i,len(t)))
                print('Source number %i, display mode %s, operation mode %s' % (I.c,display_mode,I.mode))
                print('Number keys -- select new source number')
                print('(d)ump, (m)ark components (default), mark an (o)ptical ID, (f)lag,\n    display (i)-band, display (W)ISE, go to (n)ext, (r)emove opt ID, \n    automatically separate into (t)wo or (s)ave and continue?',end=' ')
                command=input()
                if command=='d':
                    ga.dump()
                elif command=='s':
                    ga.write(name)
                    stop=True
                elif command=='f':
                    ga.flagged=True
                    ga.write(name)
                    stop=True
                elif command>='0' and command<='9':
                    I.c=int(command)
                    I.set_mode('O')
                    I.redraw()
                elif command=='n':
                    stop=True
                elif command=='i':
                    display_mode='optical'
                    break
                elif command=='r':
                    del(ga.optids[I.c])
                    ga.unchanged=False
                    I.redraw()
                elif command=='W':
                    display_mode='wise'
                    break
                elif command in ['m','o']:
                    I.set_mode(command)
                elif command=='t':
                    # Automatically select the other component of two
                    if len(mlt)!=2:
                        print('There are more than two possibilities!')
                    else:
                        cra,cdec=ga.optpos(I.c)
                        if cra is None or cdec is None:
                            print('No existing optical ID!')
                        else:
                            sep=separation(cra,cdec,mlt['ra'],mlt['dec'])
                            index=np.argmax(sep)
                            ra=mlt[index]['ra']
                            dec=mlt[index]['dec']
                            I.c=2
                            I.set_mode('O')
                            I.clickon(ra,dec)

            plt.close()

        
