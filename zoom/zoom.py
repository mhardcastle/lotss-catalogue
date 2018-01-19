#!/usr/bin/python

from astropy.table import Table,vstack
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import sys
import os
import glob
from subim import extract_subim
import matplotlib.pyplot as plt
from source_handler import Source,parsefile
from separation import separation

scale=3600.0 # scaling factor for region sizes

def find_bbox(t):
    # given a table t find the bounding box of the ellipses for the regions

    boxes=[]
    for r in t:
        a=r['Maj']/scale
        b=r['Min']/scale
        th=(r['PA']+90)*np.pi/180.0
        dx=np.sqrt((a*np.cos(th))**2.0+(b*np.sin(th))**2.0)
        dy=np.sqrt((a*np.sin(th))**2.0+(b*np.cos(th))**2.0)
        boxes.append([r['RA']-dx/np.cos(r['DEC']*np.pi/180.0),
                      r['RA']+dx/np.cos(r['DEC']*np.pi/180.0),
                      r['DEC']-dy, r['DEC']+dy])

    boxes=np.array(boxes)

    minra=np.nanmin(boxes[:,0])
    maxra=np.nanmax(boxes[:,1])
    mindec=np.nanmin(boxes[:,2])
    maxdec=np.nanmax(boxes[:,3])
    
    ra=np.mean((minra,maxra))
    dec=np.mean((mindec,maxdec))
    size=1.2*3600.0*np.max((maxdec-mindec,(maxra-minra)*np.cos(dec*np.pi/180.0)))
    return ra,dec,size

def get_mosaic_name(name):
    globst=os.environ['IMAGEDIR']+'/mosaics/'+name.rstrip()+'*'
    g=glob.glob(globst)
    if len(g)==1:
        return g[0]
    elif len(g)==0:
        raise RuntimeError('No mosaic called '+name)
    else:
        raise RuntimeError('Mosaic name ambiguous')

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
        if not(np.isnan(self.optra)):
            self.f.show_markers(self.optra,self.optdec,marker='x',facecolor='magenta',edgecolor='magenta',linewidth=3,s=1500,zorder=300,layer='ID_marker')
        else:
            try:
                self.f.remove_layer('ID_marker')
            except:
                pass
        c=[]
        for r in self.ots:
            if r['Source_Name'] in self.components:
                c.append('green')
                #self.f.show_ellipses(r['RA'],r['DEC'],r['Maj']*2/scale,r['Min']*2/scale,angle=90+r['PA'],edgecolor='green',linewidth=3,zorder=101)
            elif r['Source_Name'] in notcomponents:
                c.append('cyan')
                #self.f.show_ellipses(r['RA'],r['DEC'],r['Maj']*2/scale,r['Min']*2/scale,angle=90+r['PA'],edgecolor='cyan',linewidth=3,zorder=101)
            else:
                #self.f.show_ellipses(r['RA'],r['DEC'],r['Maj']*2/scale,r['Min']*2/scale,angle=90+r['PA'],edgecolor='red',linewidth=3,zorder=101)
                c.append('red')
        self.f.show_ellipses(ots['RA'],ots['DEC'],ots['Maj']*2/scale,ots['Min']*2/scale,angle=90+ots['PA'],edgecolor=c,linewidth=3,zorder=101,layer='Component_ellipse') 
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
                print 'Optical ID at',ra,dec
                self.redraw()
            elif self.mode=='m':
                sep=separation(ra,dec,ots['RA'],ots['DEC'])
                index=np.argmin(sep)
                name=ots[index]['Source_Name']
                if ots[index]['Source_Name'] in self.components:
                    self.components.remove(name)
                    print 'removed component',name
                else:
                    self.components.append(name)
                    print 'added component',name
                self.redraw()
            elif self.mode=='z':
                if not(np.isnan(self.oldra)):
                    self.size=separation(ra,dec,self.oldra,self.olddec)*3600
                    print 'Size measured as',self.size,'arcsec'
                self.oldra=ra
                self.olddec=dec
            else:
                raise NotImplementedError('Mode not recognised')

    def set_mode(self,mode):
        self.mode=mode
        if self.mode=='m':
            print 'Selected mark/unmark source mode'
        elif self.mode=='o':
            print 'Selected optical ID mode'
            self.optra=np.nan
            self.optdec=np.nan
        elif self.mode=='z':
            print 'Selected size mode'
            self.oldra=np.nan
            self.olddec=np.nan
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

                
if __name__=='__main__':

    imagedir=os.environ['IMAGEDIR']
    tname=sys.argv[1]
    lname=sys.argv[1].replace('.fits','-list.txt')
    t=Table.read(tname)
    # all LOFAR sources
    ot=Table.read('/data/lofar/mjh/hetdex_v4/lgzmatch/LOFAR_HBA_T1_DR1_catalog_v0.9.srl.fixed.fits')
    # large source table
    lt=ot[(ot['Total_flux']>8) & (ot['Maj']>6)]
    # galaxies if needed
    gals=Table.read(imagedir+'/wise/allwise_HETDEX_full_radec.fits')
    
    # read lists
    lines=[l.rstrip().split() for l in open(lname).readlines()]
    names=[l[0] for l in lines]
    lofarmaps=[l[1] for l in lines]
    psmaps=[l[2] for l in lines]
    wisemaps=[l[3] for l in lines]
    firstmaps=[l[4] for l in lines]

    if len(sys.argv)>1:
        clines=open('../lgz_v1/lgz_components.txt').readlines()
        ss=Source()
        for c in clines:
            bits=c.split()
            ss.add(bits[1],bits[0].rstrip())                 
    else:
        clines=open('components.txt').readlines()
        ss=Source()
        for c in clines:
            bits=c.split()
            ss.add(bits[0],bits[1].rstrip())                 

    for r in t:
        #ss.add(r['Source_Name'],r['Source_Name'])
        try:
            ss.set_opt(r['Source_Name'],r['optRA'],r['optDec'])
        except:
            try:
                ss.set_opt(r['Source_Name'],r['ID_ra'],r['ID_dec'])
            except:
                ss.set_opt(r['Source_Name'],r['LR_ra'],r['LR_dec'])
        try:
            ss.set_size(r['Source_Name'],r['Size'])
        except:
            try:
                ss.set_size(r['Source_Name'],r['LGZ_Size'])
            except:
                pass

    # now parse existing files to the structure
    dir='/data/lofar/mjh/hetdex_v4/zoom/'
    g=glob.glob(dir+'ILT*.txt')
    for f in g:
        if 'table-list' in f:
            continue
        n=f.split('/')[-1].replace('.txt','')
        print 'Parsing file for source',n
        parsefile(n,ss,dir=dir)
        ss.mdict[n]=2

    
    for i in range(len(t)):
        sourcename=names[i]
        #if os.path.isfile(sourcename+'.txt'):
        #    print sourcename,'already has a zoom file'
        #    continue
        if len(ss.get_comps(sourcename))==0:
            print sourcename,'has no components!'
            for source in ss.cdict:
                if sourcename in ss.cdict[source]:
                    break # it's already in some other source
            else:
                ss.add(sourcename,sourcename)
            if sourcename in ss.cdict[source]:
                continue
        r=t[i]
        print i,r
        assert(sourcename==r['Source_Name'])
        
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
        if os.path.isdir(lofarfile):
            lofarfile+='/mosaic.fits'

        from overlay import show_overlay

        ra,dec=r['RA'],r['DEC']

        try:
            marker_ra=r['ra']
            marker_dec=r['dec']
        except:
            marker_ra=None
            marker_dec=None

        title=sourcename

        components=ss.get_comps(sourcename)
        print 'Sourcename is',sourcename,'components is',components
        mask=np.array([False]*len(ot))
        for c in components:
            mask|=(ot['Source_Name']==c)
        ctable=ot[mask]
        print ctable
        # resize the image to look for interesting neighbours
        iter=0
        dist=180
        while True:
            startra,startdec=ra,dec
            tcopy=lt
            tcopy['dist']=np.sqrt((np.cos(dec*np.pi/180.0)*(tcopy['RA']-ra))**2.0+(tcopy['DEC']-dec)**2.0)*3600.0
            tcopy=tcopy[tcopy['dist']<dist]
            if len(tcopy)==0:
                dist*=1.2
                continue
            
            # make sure the original source is in there
            tcopy=vstack((ctable,r))

            ra=np.mean(tcopy['RA'])
            dec=np.mean(tcopy['DEC'])

            if startra==ra and startdec==dec:
                break
            iter+=1
            if iter==10:
                break

        # now find the bounding box of the resulting collection
        ra,dec,size=find_bbox(tcopy)
        print tcopy,size

        if np.isnan(size):
            ra=r['RA']
            dec=r['DEC']
            size=60

        if size>300:
            size=300.0
        if size<60:
            size=60.0
        size=(int(0.5+size/10))*10
        # zoom out
        size*=2

        size/=3600.0

        '''
        gals=Table.read(imagedir+'/tier1_ps1_wise_hetdex.fits')
        gals['raMean'].name='ra'
        gals['decMean'].name='dec'
        pg=gals[(np.abs(gals['ra']-ra)<size) & (np.abs(gals['dec']-dec)<size)]
        del(gals)

        '''
        pwg=gals[(np.abs(gals['ra']-ra)<size) & (np.abs(gals['dec']-dec)<size)]

        ots=ot[separation(ra,dec,ot['RA'],ot['DEC'])<(size*2)]

        #pshdu=extract_subim(imagedir+'/downloads/'+psmaps[i],ra,dec,size*2,hduid=1)
        print 'Lofarfile is',lofarfile
        lhdu=extract_subim(lofarfile,ra,dec,size*2)
        firsthdu=extract_subim(imagedir+'/downloads/'+firstmaps[i],ra,dec,size*2)
        whdu=extract_subim(imagedir+'/downloads/'+wisemaps[i],ra,dec,size)
        try:
            peak==r['Peak_flux']/1000.0
        except:
            peak=None
        
        f=show_overlay(lhdu,whdu,ra,dec,size,firsthdu=firsthdu,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,block=False,interactive=False,drlimit=8000,peak=peak,plotpos=[(pwg,'+')])

        ora,odec=ss.odict[sourcename]
        components=ss.get_comps(sourcename)
        notcomponents=ss.get_ncomps(sourcename)
        print 'components is',components
        I=Interactive(f,ora,odec,ots,components,notcomponents)
        fig=plt.gcf()
        fig.canvas.mpl_connect('button_press_event', I.onclick)

        plt.ion()
        plt.show(block=False)
        
        stop=False
        while not(stop):
            print '(d)rop source, (m)ark components (default), mark an (o)ptical ID,\n   mark a si(z)e, set (b)lend, go to (n)ext or (s)ave and continue?',
            command=raw_input()
            if command=='s':
                stop=True
                I.write(r['Source_Name'])
            elif command=='d':
                stop=True
                I.delete(r['Source_Name'])
            elif command=='n':
                stop=True
            elif command in ['m','o','z']:
                I.set_mode(command)
            elif command=='b':
                I.blend=True
            elif command=='p':
                continue
            else:
                print 'Command not recognised!'

        if command=='s':
            ss.set_components(sourcename,I.components)
            ss.set_opt(sourcename,I.optra,I.optdec)
            ss.set_size(sourcename,I.size)
        elif command=='d':
            ss.delete_source(sourcename)

        # options are:
        # drop source entirely
        # mark/unmark components
        # mark a size (for if components get this completely wrong)?
        # mark an optical ID (or just unmark existing if nothing selected)

        plt.close()
