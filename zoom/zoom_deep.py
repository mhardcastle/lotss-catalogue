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
from make_catalogue import Source,make_structure,generate_table
from separation import separation
from image_utils import find_bbox,get_mosaic_name
from get_fits import save_fits
from overlay import show_overlay
scale=1.0 # units of catalogue are degrees

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
                print 'Optical ID at',ra,dec
                self.redraw()
            elif self.mode=='m':
                sep=separation(ra,dec,self.ots['RA'],self.ots['DEC'])
                index=np.argmin(sep)
                name=self.ots[index]['Source_Name']
                if self.ots[index]['Source_Name'] in self.components:
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
            elif self.mode=='i':
                sep=separation(ra,dec,self.ots['RA'],self.ots['DEC'])
                index=np.argmin(sep)
                print self.ots[index]
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
            self.redraw()
        elif self.mode=='z':
            print 'Selected size mode'
            self.oldra=np.nan
            self.olddec=np.nan
        elif self.mode=='i':
            print 'Selected component inspector mode'
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

    dir=os.getcwd()
    field=os.path.basename(dir)
    print 'field is',field
    if not os.path.isdir('zoom'):
        os.mkdir('zoom') # store text files in here

    s=make_structure(field,warn=True)
    print 'Reading files'
    if field=='bootes':

        sourcecat = "/beegfs/lofar/deepfields/Bootes_LR/new_fdeep_matches/Bootes_ML_RUN_fin_overlap_srl_workflow_th.fits"
        optcat='/beegfs/lofar/deepfields/Bootes_merged_optical/Bootes_MASTER_opt_spitzer_merged.fits'
        opt=Table.read(optcat)
        flt=(opt['FLAG_DEEP']>0)
        flt&=(opt['FLAG_OVERLAP']>0)
        opt=opt[flt]
        idname='ID'
        lofarfile=fits.open('/beegfs/lofar/deepfields/Bootes_LOFAR/image_full_ampphase_di_m.NS_shift.int.facetRestored.blanked.scaled.fits')
        spitzerfile=fits.open('/beegfs/lofar/deepfields/Bootes_optical/SDWFS/I2_bootes.v32.fits')
        ibandfile=fits.open('/beegfs/lofar/deepfields/Bootes_merged_optical/Bootes_iband.fits')

    elif field=='lockman':
        sourcecat = "/beegfs/lofar/deepfields/Lockman_LR/updated_LR_cols/LH_ML_RUN_fin_overlap_srl_workflow_th.fits"
        optcat='/beegfs/lofar/deepfields/Lockman_edited_cats/optical/LH_MASTER_opt_spitzer_merged_forLGZ.fits'
        opt=Table.read(optcat)
        flt=(opt['FLAG_OVERLAP']==3)
        opt=opt[flt]
        idname='NUMBER'
        lofarfile=fits.open('/beegfs/lofar/deepfields/Lockman_LOFAR/image_full_ampphase_di_m.NS_shift.int.facetRestored.blanked.scaled.fits')
        spitzerfile=fits.open('/beegfs/lofar/deepfields/Lockman/LH_4d5band_old.fits')
        rbandfile=fits.open('/beegfs/lofar/deepfields/Lockman/LH_rband.fits')

    elif field=='en1':
        sourcecat = "/beegfs/lofar/deepfields/ELAIS_N1_LR/new_optcat_matches/EN1_ML_RUN_fin_overlap_srl_workflow_th.fits"
        optcat='/beegfs/lofar/deepfields/ELAIS_N1_optical/catalogues/correct_merging/add_uncat/EN1_merged_pos.fits'
        opt=Table.read(optcat)
        flt=(opt['FLAG_OVERLAP']==7)
        opt=opt[flt]
        idname='NUMBER'
        lofarfile=fits.open('/beegfs/lofar/deepfields/ELAIS-N1_LOFAR/image_full_ampphase_di_m.NS_shift.int.facetRestored.fits')
        spitzerfile=fits.open('/beegfs/lofar/deepfields/ELAIS_N1_optical/optical_images/sw2band/EL_EN1_sw2band.fits')
        ibandfile=fits.open('/beegfs/lofar/deepfields/ELAIS_N1_optical/optical_images/iband/EL_EN1_iband.fits')

    gals=opt
    if 'ALPHA_J2000' in gals.colnames:
        gals['ALPHA_J2000'].name='ra'
    if 'DELTA_J2000' in gals.colnames:
        gals['DELTA_J2000'].name='dec'

    # all LOFAR sources
    # we generate this on the fly from the s structure, since otherwise it doesn't take any deblending etc into account
    columns=[('Source_Name',None),('RA',None),('DEC',None),('E_RA',None),('E_DEC',None),('Total_flux',None),('E_Total_flux',None),('Peak_flux',None),('E_Peak_flux',None),('S_Code',None),('Maj',np.nan),('Min',np.nan),('PA',np.nan),('E_Maj',np.nan),('E_Min',np.nan),('E_PA',np.nan),('DC_Maj',np.nan),('DC_Min',np.nan),('DC_PA',np.nan),('Created',None),('Parent',None)]
    ot=generate_table(s.cd,columns,keep_deleted=True)
    # large source table
    lt=ot[(ot['Total_flux']>0.005) & (ot['Maj']>10/3600.0)]

    os.chdir('zoom')
    sourcelist=[]
    try:
        sourcename=sys.argv[1]
    except:
        sourcename=None
    if sourcename is None:
        for sourcename in s.sd:
            if 'Deleted' in s.sd[sourcename]:
                continue
            if ( ('Zoom_prob' in s.sd[sourcename] and s.sd[sourcename]['Zoom_prob']>0.5) or
                 ('Imagemissing_prob' in s.sd[sourcename] and s.sd[sourcename]['Imagemissing_prob']>0.5) or
                 ('Hostbroken_prob' in s.sd[sourcename] and s.sd[sourcename]['Hostbroken_prob']>0.5) ):
            
                sourcelist.append(sourcename)
    else:
        if sourcename.startswith('ILTJ'):
            sourcelist=[sourcename]
        else:
            # maybe this is a file with a list of sources
            if not os.path.isfile(sourcename):
                raise RuntimeError('Cannot parse entity on command line')
            else:
                sourcelist=[l.rstrip() for l in open(sourcename).readlines()]

    sourcelist+=s.zoomneeded
    mode='spitzer'
    for sourcename in sourcelist:
        if len(sourcelist)>1 and os.path.isfile(sourcename+'.txt'):
            print sourcename,'already has a zoom file'
            continue
        if len(s.get_comps(sourcename))==0:
            print sourcename,'has no components!'
            continue
            
        ra,dec=s.sd[sourcename]['RA'],s.sd[sourcename]['DEC']

        marker_ra=None
        marker_dec=None

        title=sourcename

        components=s.get_comps(sourcename)
        print 'Sourcename is',sourcename,'components is',components
        mask=np.array([False]*len(ot))
        for c in components:
            mask|=(ot['Source_Name']==c)
        ctable=ot[mask]
        print ctable
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
        print tcopy,size

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
        initfactor=2.0
        scalefactor=initfactor
        overlaygals=True
        while True:
            print 'Scale factor is',scalefactor
            if overlaygals:
                pwg=gals[(np.abs(gals['ra']-ra)<size*scalefactor/np.cos(dec*np.pi/180.0)) & (np.abs(gals['dec']-dec)<size*3)]
            else:
                pwg=None

            ots=ot[separation(ra,dec,ot['RA'],ot['DEC'])<(size*scalefactor)]

            #pshdu=extract_subim(imagedir+'/downloads/'+psmaps[i],ra,dec,size*2,hduid=1)
            print 'Lofarfile is',lofarfile
            lhdu=extract_subim(lofarfile,ra,dec,size*scalefactor)
            if mode=='spitzer':
                whdu=extract_subim(spitzerfile,ra,dec,size*scalefactor)
            else:
                whdu=extract_subim(ibandfile,ra,dec,size*scalefactor)
            try:
                peak==r['Peak_flux']/1000.0
            except:
                peak=None

            try:
                f=show_overlay(lhdu,whdu,ra,dec,size*scalefactor/initfactor,coords_color='red',coords_ra=s.sd[sourcename]['RA'],coords_dec=s.sd[sourcename]['DEC'],coords_lw=3,lw=2,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,block=False,interactive=False,drlimit=8000,peak=peak,plotpos=[(pwg,'+')])
            except Exception as e:
                print '*** Error making image ***'
                print e
                raw_input()
                stop=True
                break

            ora,odec=s.optid(sourcename)

            components=s.get_comps(sourcename)
            notcomponents=s.get_ncomps(sourcename) # list of all components in other sources
            print 'components is',components
            I=Interactive(f,ora,odec,ots,components,notcomponents)
            fig=plt.gcf()
            fig.canvas.mpl_connect('button_press_event', I.onclick)

            plt.ion()
            plt.show(block=False)

            stop=False
            while not(stop):
                print '(d)rop source, (m)ark components (default), mark an (o)ptical ID,\n   mark a si(z)e, set (b)lend, save (f)its, go to (n)ext, (Z)oom out, \n   (i)nspect, (T)oggle galaxy overlay, change opt (I)mage or (s)ave and continue?',
                command=raw_input()
                if command=='s':
                    stop=True
                    I.write(sourcename)
                elif command=='d':
                    stop=True
                    I.delete(sourcename)
                elif command=='n':
                    stop=True
                elif command in ['m','o','z','i']:
                    I.set_mode(command)
                elif command=='b':
                    I.blend=True
                elif command=='p':
                    continue
                elif command=='f':
                    lhdu.writeto(sourcename+'.fits',overwrite=True)
                    print 'Saved as',sourcename+'.fits'
                    os.system('ds9 '+sourcename+'.fits &')
                elif len(command)>0 and command[0]=='Z':
                    try:
                        scalefactor=int(command[1:])
                    except:
                        scalefactor*=2
                    plt.close()
                    break
                elif command=='T':
                    overlaygals=not(overlaygals)
                    plt.close()
                    break
                elif command=='I':
                    if mode=='spitzer':
                        mode='optical'
                    else:
                        mode='spitzer'
                    plt.close()
                    break
                else:
                    print 'Command not recognised!'
            if command=='s':
                s.set_components(sourcename,I.components)
                s.set_opt(sourcename,I.optra,I.optdec)
                s.set_size(sourcename,I.size)
            elif command=='d':
                s.delete_source(sourcename,'Marked as LGZ artefact')
            
            if stop:
                break # out of outer while for scalefactor reset
            
        # options are:
        # drop source entirely
        # mark/unmark components
        # mark a size (for if components get this completely wrong)?
        # mark an optical ID (or just unmark existing if nothing selected)

        plt.close()
