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

scale=3600.0

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
        outfile=open(name+'.txt','w')
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
        print 'Components:'
        for i,cc in enumerate(self.components):
            print i,cc,self.assocs[cc],self.ots[i]['Total_flux']
        print 'Optical IDs'
        for k in self.optids:
            print k,self.optids[k]
        
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
            print '\nOptical ID at',ra,dec
            if self.mode=='o': self.redraw()
        if self.mode=='m' or self.mode=='O':
            sep=separation(ra,dec,self.ots['RA'],self.ots['DEC'])
            index=np.argmin(sep)
            name=self.ots[index]['Gaus_id']
            if self.ga.ismember(self.c,name):
                self.ga.remove(name)
                print '\nremoved component',name,'from source',self.c
            else:
                self.ga.add(self.c,name)
                print '\nadded component',name,'to source',self.c
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
            print 'Selected mark/unmark source mode'
        elif self.mode=='o':
            print 'Selected optical ID mode'
        elif self.mode=='O':
            print 'Selected select-and-ID mode'
        else:
            raise NotImplementedError('Mode not recognised')

if __name__=='__main__':

    imagedir=os.environ['IMAGEDIR']
    tname=sys.argv[1]
    lname=sys.argv[1].replace('.fits','-list.txt')
    t=Table.read(tname)
    # read lists
    lines=[l.rstrip().split() for l in open(lname).readlines()]
    names=[l[0] for l in lines]
    lofarmaps=[l[1] for l in lines]
    psmaps=[l[2] for l in lines]
    wisemaps=[l[3] for l in lines]
    firstmaps=[l[4] for l in lines]
    # galaxies if needed
    #gals=Table.read(imagedir+'/wise/allwise_HETDEX_full_radec.fits')
    gals=Table.read('/data/lofar/mjh/hetdex_ps1_allwise_radec.fits')
    ct=Table.read('LOFAR_HBA_T1_DR1_merge_ID_v1.1.comp.fits')
    t['Source_Name']=[s.rstrip() for s in t['Source_Name']]
    ct['Source_Name']=[s.rstrip() for s in ct['Source_Name']]
    #gt=Table.read('LOFAR_HBA_T1_DR1_catalog_v0.9.gaus.fixed.fits')
    gt=Table.read('lofar_gaus_pw.fixed.fits')
    display_mode='WISE'

    for i,r in enumerate(t):
        name=r['Source_Name']
        if os.path.isfile(name+'.txt'):
            print name,'already has a blend file'
            continue
        gc=0
        gl=[]
        print '>>%s<<' % name
        ctf=(ct['Source_Name']==name)
        ctfl=ct[ctf]
        print '... has',len(ctfl),'components'
        if len(ctfl)==0:
            continue
        for j,c in enumerate(ctfl):
            print '    Component',j,'has type',c['S_Code']
            gtf=(gt['Source_Name']==c['Component_Name'])
            gtfl=gt[gtf]
            print '    ... and contains',len(gtfl),'Gaussians'
            gl.append(gtfl)
        rgt=vstack(gl)
        print '... altogether',len(rgt),'Gaussians'
        mlt=rgt[~np.isnan(rgt['ra'])]
        try:
            mosaic=r['Mosaic_ID']
        except:
            mosaic=None
        print 'Mosaic id is',mosaic
        if mosaic is not None:
            try:
                lofarfile=get_mosaic_name(mosaic)
            except RuntimeError:
                print 'Could not get mosaic ID from',r['Mosaic_ID']
                mosaic=None
        if mosaic is None:
            lofarfile=os.environ['IMAGEDIR']+'/'+lofarmaps[i]
        title=name
        ra,dec,size=find_bbox(rgt)
        size*=1.5
        size/=3600.0
        print 'Lofarfile is',lofarfile
        lhdu=extract_subim(lofarfile,ra,dec,size*2)
        firsthdu=extract_subim(imagedir+'/downloads/'+firstmaps[i],ra,dec,size*2)
        try:
            peak==r['Peak_flux']/1000.0
        except:
            peak=None
        pwg=gals[(np.abs(gals['ra']-ra)<size) & (np.abs(gals['dec']-dec)<size)]
        
        ora=r['ID_ra']
        odec=r['ID_dec']
        ga=GaussAssoc(rgt,ora,odec)

        stop=False
        while not(stop):
            if display_mode=='WISE':
                ohdu=extract_subim(imagedir+'/downloads/'+wisemaps[i],ra,dec,size)
            else:
                ohdu=extract_subim(imagedir+'/downloads/'+psmaps[i],ra,dec,size*2,hduid=1)
            f=show_overlay(lhdu,ohdu,ra,dec,size,firsthdu=firsthdu,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,no_labels=True,title=title,block=False,interactive=False,drlimit=8000,peak=peak,plotpos=[(pwg,'+'),(mlt,'h','red')])

            I=Interactive(f,ga)
            fig=plt.gcf()
            fig.canvas.mpl_connect('button_press_event', I.onclick)

            plt.ion()
            plt.show(block=False)

            while not(stop):
                print 'Source number %i, display mode %s, operation mode %s' % (I.c,display_mode,I.mode)
                print 'Number keys -- select new source number'
                print '(d)ump, (m)ark components (default), mark an (o)ptical ID, (f)lag,\n    display (p)an-STARRS, display (w)ISE, go to (n)ext, (r)emove opt ID, \n    automatically separate into (t)wo or (s)ave and continue?',
                command=raw_input()
                if command=='d':
                    ga.dump()
                    print pwg
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
                elif command=='p':
                    display_mode='Pan-STARRS'
                    break
                elif command=='r':
                    del(ga.optids[I.c])
                    ga.unchanged=False
                    I.redraw()
                elif command=='w':
                    display_mode='WISE'
                    break
                elif command in ['m','o']:
                    I.set_mode(command)
                elif command=='t':
                    # Automatically select the other component of two
                    if len(mlt)!=2:
                        print 'There are more than two possibilities!'
                    else:
                        cra,cdec=ga.optpos(I.c)
                        if cra is None or cdec is None:
                            print 'No existing optical ID!'
                        else:
                            sep=separation(cra,cdec,mlt['ra'],mlt['dec'])
                            index=np.argmax(sep)
                            ra=mlt[index]['ra']
                            dec=mlt[index]['dec']
                            I.c=2
                            I.set_mode('O')
                            I.clickon(ra,dec)

            plt.close()

        
