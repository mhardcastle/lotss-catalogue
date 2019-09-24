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

scale=1.0 # units of catalogue are degrees

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

def preselect_read(dir,tname):
    t=Table.read(tname)
    group=[]
    source=[]
    filt=[False]*len(t)
    wc=dir+'/*/workflow.txt'
    print 'wild card is',wc
    g=glob.glob(wc)
    for f in g:
        print f
        lines=open(f).readlines()
        for l in lines:
            l=l.rstrip()
            bits=l.split(',')
            group.append(int(bits[0]))
            source.append(bits[1])
    for g,s in zip(group,source):
        print g,s
        if g==6 or g==7:
            i=np.argmax(t['Source_Name']==s)
            filt[i]=True
    print 'Analysing',np.sum(filt),'preselected objects'
    return t[filt]
        
if __name__=='__main__':

    dir=os.getcwd()
    field=os.path.basename(dir)
    print 'field is',field
    if not os.path.isdir('blend'):
        os.mkdir('blend') # store text files in here

    print 'Reading data...'

    if field=='bootes':
        sourcecat = "/beegfs/lofar/deepfields/Bootes_LR/new_fdeep_matches/Bootes_ML_RUN_fin_overlap_srl_workflow_th.fits"
        t=preselect_read('/beegfs/lofar/deepfields/Bootes_preselect',sourcecat)

        optcat='/beegfs/lofar/deepfields/Bootes_merged_optical/Bootes_MASTER_opt_spitzer_merged.fits'
        opt=Table.read(optcat)
        flt=(opt['FLAG_DEEP']>0)
        flt&=(opt['FLAG_OVERLAP']>0)
        opt=opt[flt]
        idname='ID'
        lofarfile=fits.open('/beegfs/lofar/deepfields/Bootes_LOFAR/image_full_ampphase_di_m.NS_shift.int.facetRestored.blanked.scaled.fits')
        spitzerfile=fits.open('/beegfs/lofar/deepfields/Bootes_optical/SDWFS/I2_bootes.v32.fits')
        ibandfile=fits.open('/beegfs/lofar/deepfields/Bootes_merged_optical/Bootes_iband.fits')
        st=Table.read('/beegfs/lofar/deepfields/Bootes_LR/new_fdeep_matches/Bootes_ML_RUN_fin_overlap_srl_workflow_th.fits')
        gt=Table.read('/beegfs/lofar/deepfields/Bootes_LR/new_fdeep_matches/Bootes_ML_RUN_fin_overlap_gaul_workflow_th.fits')

    elif field=='lockman':
        sourcecat = "/beegfs/lofar/deepfields/Lockman_LR/LH_ML_RUN_fin_overlap_srl_workflow_th.fits"
        t=preselect_read('/beegfs/lofar/deepfields/Lockman_preselect',sourcecat)
        optcat='/beegfs/lofar/deepfields/Lockman_edited_cats/optical/LH_MASTER_opt_spitzer_merged_forLGZ.fits'
        opt=Table.read(optcat)
        flt=(opt['FLAG_OVERLAP']==3)
        opt=opt[flt]
        idname='NUMBER'
        lofarfile=fits.open('/beegfs/lofar/deepfields/Lockman_LOFAR/image_full_ampphase_di_m.NS_shift.int.facetRestored.blanked.scaled.fits')
        spitzerfile=fits.open('/beegfs/lofar/deepfields/Lockman/LH_4d5band.fits')
        rbandfile=fits.open('/beegfs/lofar/deepfields/Lockman/LH_rband.fits')
        st=Table.read('/beegfs/lofar/deepfields/Lockman_LR/updated_LR_cols/LH_ML_RUN_fin_overlap_srl_workflow_th.fits')
        gt=Table.read('/beegfs/lofar/deepfields/Lockman_LR/updated_LR_cols/LH_ML_RUN_fin_overlap_gaul_workflow_th.fits')

    elif field=='en1':
        sourcecat = "/beegfs/lofar/deepfields/ELAIS_N1_LR/new_optcat_matches/EN1_ML_RUN_fin_overlap_srl_workflow_th.fits"
        t=preselect_read('/beegfs/lofar/deepfields/ELAIS-N1_preselect',sourcecat)
        optcat='/beegfs/lofar/deepfields/ELAIS_N1_optical/catalogues/correct_merging/EN1_MASTER_opt_spitzer_merged_cedit_apcorr.fits'
        opt=Table.read(optcat)
        flt=(opt['FLAG_OVERLAP']==7)
        opt=opt[flt]
        idname='NUMBER'
        lofarfile=fits.open('/beegfs/lofar/deepfields/ELAIS-N1_LOFAR/image_full_ampphase_di_m.NS_shift.int.facetRestored.fits')
        spitzerfile=fits.open('/beegfs/lofar/deepfields/ELAIS_N1_optical/optical_images/sw2band/EL_EN1_sw2band.fits')
        ibandfile=fits.open('/beegfs/lofar/deepfields/ELAIS_N1_optical/optical_images/iband/EL_EN1_iband.fits')
        st=Table.read('/beegfs/lofar/deepfields/ELAIS_N1_LR/new_optcat_matches/EN1_ML_RUN_fin_overlap_srl_workflow_th.fits')
        gt=Table.read('/beegfs/lofar/deepfields/ELAIS_N1_LR/new_optcat_matches/EN1_ML_RUN_fin_overlap_gaul_workflow_th.fits')
    else:
        print 'Not in correct working directory'
        sys.exit(1)
        
    gals=opt
    if 'ALPHA_J2000' in gals.colnames:
        gals['ALPHA_J2000'].name='ra'
    if 'DELTA_J2000' in gals.colnames:
        gals['DELTA_J2000'].name='dec'
        
    display_mode='Spitzer'
    
    for i,r in enumerate(t):
        name=r['Source_Name']
        print 'Source is',name
        # work out whether there are sources selected in LGZ
            
        sts=st[st['Source_Name']==name]
        if len(sts)==0:
            print name,'does not exist in source table, skipping'
            continue
        r2=sts[0] # corresponding row in source table
        if os.path.isfile('blend/'+name+'.txt'):
            print name,'already has a blend file'
            continue
        id=r2['Source_id'] # unique in this field!
        rgt=gt[gt['Source_id']==id]
        print 'Contains',len(rgt),'Gaussians'
        mlt=rgt[~np.isnan(rgt['lr_ra_fin'])]
        mlt['ra']=mlt['lr_ra_fin']
        mlt['dec']=mlt['lr_dec_fin']
        title=name+' (%s)' % display_mode
        ra,dec,size=find_bbox(rgt,scale=scale)
        size*=1.5
        size/=3600.0
        print 'Lofarfile is',lofarfile
        lhdu=extract_subim(lofarfile,ra,dec,size*2)
        try:
            peak==r2['Peak_flux']/1000.0
        except:
            peak=None
        pwg=gals[(np.abs(gals['ra']-ra)<size) & (np.abs(gals['dec']-dec)<size)]

        ora=r2['lr_ra_fin']
        odec=r2['lr_dec_fin']
        ga=GaussAssoc(rgt,ora,odec)

        stop=False
        while not(stop):
            if display_mode=='Spitzer':
                ohdu=extract_subim(spitzerfile,ra,dec,size)
            else:
                ohdu=extract_subim(ibandfile,ra,dec,size)
            f=show_overlay(lhdu,ohdu,ra,dec,size,coords_color='red',coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,no_labels=True,title=title,block=False,interactive=False,drlimit=8000,peak=peak,plotpos=[(pwg,'+'),(mlt,'h','red')])

            I=Interactive(f,ga)
            fig=plt.gcf()
            fig.canvas.mpl_connect('button_press_event', I.onclick)

            plt.ion()
            plt.show(block=False)

            while not(stop):
                print 'Source number %i, display mode %s, operation mode %s' % (I.c,display_mode,I.mode)
                print 'Number keys -- select new source number'
                print '(d)ump, (m)ark components (default), mark an (o)ptical ID, (f)lag,\n    display (i)-band, display (S)pitzer, go to (n)ext, (r)emove opt ID, \n    automatically separate into (t)wo or (s)ave and continue?',
                command=raw_input()
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
                elif command=='S':
                    display_mode='Spitzer'
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

        
