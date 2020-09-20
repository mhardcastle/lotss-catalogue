#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
from astropy.table import Table,vstack
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import sys
import os
import glob
from separation import separation
from subim import extract_subim
from image_utils import find_bbox,get_mosaic_name
from overlay import show_overlay
from shapely.geometry import Polygon

def ellipse_poly(x0,y0,a,b,pa,n=200):
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

def ellipse(r,ra,dec):
    return ellipse_poly((ra-r['RA'])*np.cos(np.pi*dec/180.0)*3600,(r['DEC']-dec)*3600,r['Maj'],r['Min'],r['PA'])

def intersects(r1,r2):
    return r1['ellipse'].intersects(r2['ellipse'])


scale=3600.0 # scaling factor for region sizes

if __name__=='__main__':

    logfile=open('log.txt','w')
    imagedir=os.environ['IMAGEDIR']
    tname=sys.argv[1]
    lname=sys.argv[1].replace('.fits','-list.txt')
    t=Table.read(tname)
    # Annotated PyBDSF table
    ot=Table.read(os.environ['LOTSS_COMPONENT_CATALOGUE'])
    if 'Component_Name' in ot.colnames:
        cname='Component_Name'
    else:
        cname='Source_Name'
    
    # large source table for neighbours
    #lt=ot[(ot['Total_flux']>3) & (ot['Maj']>8)]
    filt=ot['Maj']>8
    filt&=ot['Total_flux']>1
    filt&=ot['Peak_flux']>2*ot['Isl_rms']
    lt=ot[filt]
    
    #print lt
    # read lists
    lines=[l.rstrip().split() for l in open(lname).readlines()]
    lofarmaps=[l[1] for l in lines]
    psmaps=[l[-1] for l in lines if len(l)>1]
    #wisemaps=[l[3] for l in lines]
    #firstmaps=[l[4] for l in lines]

    start=int(sys.argv[2])
    try:
        end=int(sys.argv[3])+1
        if end>len(t):
            end=len(t)
    except:
        end=start+1
    for i in range(start,end):
        r=t[i]
        #print r
        try:
            sourcename=r['Source_Name']
        except KeyError:
            # try to fix up old-format files
            t['Source_id'].name='Source_Name'
            sourcename=r['Source_Name']
        sourcename=sourcename.rstrip()

        csimage=sourcename+'_CS.png'
        psimage=sourcename+'_PS.png'
        pspimage=sourcename+'_PSp.png'
        manifestname=sourcename+'-manifest.txt'
        #print lofarmaps[i],psmaps[i]

        try:
            mosaic=r['Mosaic_ID']
        except:
            mosaic=None
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
            if not os.path.isfile(lofarfile):
                lofarfile=lofarfile.replace('mosaic.fits','mosaic-blanked.fits')
            

        #print lofarfile
        if '.fits' not in lofarfile:
            raise RuntimeError('Could not find mosaic fits file')
            #print 'Could not find mosaic fits file, skipping {}'.format(r['Source_Name'])
            #continue 

        if os.path.isfile(manifestname):
            print 'Selected output file exists already'
            continue

        ra,dec=r['RA'],r['DEC']
        print 'Original RA DEC is',ra,dec
        
        #try:
        #    marker_ra=r['ra']
        #    marker_dec=r['dec']
        #except:
        marker_ra=None
        marker_dec=None

        title=None

        # new part of the algorithm: first look for any ellipses that
        # intersect the current one(s) (using ot)

        iter=0
        checktable=Table(r)
        checktable['ellipse']=None
        checktable['ellipse'][0]=ellipse(r,ra,dec)
        ctlen=1
        while True:
            tcopy=ot
            tcopy['dist']=np.sqrt((np.cos(dec*np.pi/180.0)*(tcopy['RA']-ra))**2.0+(tcopy['DEC']-dec)**2.0)*3600.0
            filt=tcopy['dist']<180
            filt&=tcopy['dist']>0
            tcopy=tcopy[filt]
            if len(tcopy)==0:
                print 'No neighbours found!'
                break
            tcopy['ellipse']=None
            for j in range(len(tcopy)):
                tcopy[j]['ellipse']=ellipse(tcopy[j],ra,dec)
            intersect=[]
            for row in checktable:
                filt=[]
                for row2 in tcopy:
                    if row[cname]==row2[cname]:
                        result=False
                    else:
                        result=intersects(row,row2)
                    print 'Check intercept between',row[cname],'and',row2[cname],'result is',result
                    
                    filt.append(result)
            print 'Filt is',filt
            itable=tcopy[filt]
            print 'Length of checktable is',len(checktable)
            print 'Length of itable is',len(itable)
            interim=vstack((checktable,itable))
            print 'Length of interim table is',len(interim)
            names=list(set(interim[cname]))
            print 'Names are',names
            filt=[False]*len(interim)
            for n in names:
                pos=np.argmax(interim[cname]==n)
                filt[pos]=True
            checktable=interim[filt]
            print 'New checktable length is',len(checktable)
            if len(checktable)==ctlen:
                print 'Checktable length converged!'
                break
            ctlen=len(checktable)
            iter+=1
            if iter==10:
                print 'Not converged!'
                break

        print checktable
        '''
        import matplotlib.pyplot as plt
        for r in checktable:
            x,y=r['ellipse'].exterior.xy
            plt.plot(x,y)
        plt.savefig('poly.pdf')
        '''
        #ra,dec,size=find_bbox(checktable)
        
        # resize the image to look for interesting neighbours (using lt)
        iter=0
        flux=r['Total_flux']
        while True:
            startra,startdec=ra,dec
            tcopy=lt
            tcopy['dist']=np.sqrt((np.cos(dec*np.pi/180.0)*(tcopy['RA']-ra))**2.0+(tcopy['DEC']-dec)**2.0)*3600.0
            filt=tcopy['dist']<180
            filt&=tcopy['Total_flux']>flux/3.0 # look for similar flux
            tcopy=tcopy[filt]
            print 'Iter',iter,'found',len(tcopy),'neighbours'

            # include checktable
            interim=vstack((tcopy,checktable))
            # de-dupe again
            names=list(set(interim[cname]))
            print 'Names are',names
            filt=[False]*len(interim)
            for n in names:
                pos=np.argmax(interim[cname]==n)
                filt[pos]=True
            tcopy=interim[filt]
            
            ra=np.mean(tcopy['RA'])
            dec=np.mean(tcopy['DEC'])
            
            newra,newdec,size=find_bbox(tcopy)
            print '     size is',size

            if startra==ra and startdec==dec:
                break
            iter+=1
            if iter==10:
                break

        # now find the bounding box of the resulting collection
        ra,dec,size=find_bbox(tcopy)

        if np.isnan(size):
            ra=r['RA']
            dec=r['DEC']
            size=60

        if size>300.0:
            # revert just to original
            ra,dec=r['RA'],r['DEC']
            tcopy=checktable
            ra,dec,size=find_bbox(tcopy)

        if size>300:
            size=300.0
        if size<60:
            size=60.0
        size=(int(0.5+size/10))*10
        print 'final size is',size
        print 'final RA DEC is',ra,dec
        
        size/=3600.0

        seps=separation(ra,dec,ot['RA'],ot['DEC'])
        ots=ot[seps<size*2]

        #ots=ots[ots['Source_Name']!=""] # removes artefacts
        #ots.write(sourcename+'-ellipses.fits',overwrite=True)
        ls=[]
        for nr in ots:
            print nr[cname],sourcename
            if nr[cname].rstrip()==sourcename:
                ls.append('solid')
            else:
                ls.append('dashed')
        print ls

        if 'Isl_rms' in r.colnames:
            rms=r['Isl_rms']/1000.0
        else:
            rms=None

        optfile=imagedir+'/downloads/'+psmaps[i]
        print 'Using optical image',optfile
        pshdu=fits.open(optfile)
        if pshdu[0].header['NAXIS']==0:
            print '*** No optical image! ***'
            logfile.write('*** No optical %s %f %f ***\n' % (sourcename,r['RA'],r['DEC']))
            continue
        # nan-blank
        pshdu[0].data=np.where(pshdu[0].data>8,np.nan,pshdu[0].data)
        
        lhdu=extract_subim(lofarfile,ra,dec,size/1.5)
        lhdu.writeto('lofar.fits',clobber=True)
        #firsthdu=extract_subim(imagedir+'/downloads/'+firstmaps[i],ra,dec,size)
        try:
            peak==r['Peak_flux']/1000.0
        except:
            peak=None

        show_overlay(lhdu,lhdu,ra,dec,size,firsthdu=None,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_lw=3,lw=1,save_name=csimage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,peak=peak,plot_coords=False,show_grid=False,lw_ellipse=3,ellipse_style=ls,ellipse_color='cyan',noisethresh=1.5,drlimit=1000,rms_use=rms,lofarlevel=2.5,logfile=logfile,sourcename=sourcename,vmax_cap=0.5,show_lofar=False,cmap='inferno',lofar_colorscale=True)

        show_overlay(lhdu,pshdu,ra,dec,size,firsthdu=None,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_lw=3,lw=1,save_name=psimage,no_labels=True,marker_ra=marker_ra,marker_dec=marker_dec,marker_lw=3,marker_color='cyan',title=title,peak=peak,plot_coords=False,show_grid=False,lw_ellipse=3,ellipse_style=ls,ellipse_color='cyan',noisethresh=1.5,drlimit=1000,rms_use=rms,lofarlevel=2.5,logfile=logfile,sourcename=sourcename,vmax_cap=0.5)
        
        show_overlay(lhdu,pshdu,ra,dec,size,overlay_cat=ots,overlay_scale=scale,coords_color='red',coords_lw=3,lw=2,show_lofar=False,save_name=pspimage,no_labels=True,title=title,peak=peak,plot_coords=False,show_grid=False,lw_ellipse=3,ellipse_style=ls,ellipse_color='cyan',noisethresh=1.5,drlimit=1000,rms_use=rms,lofarlevel=2.5,vmax_cap=0.5)

        #wiseimage=sourcename+'_W.png'
        #firsthdu=extract_subim(imagedir+'/downloads/'+firstmaps[i],ra,dec,size)
        #whdu=extract_subim(imagedir+'/downloads/'+wisemaps[i],ra,dec,size)
        #show_overlay(lhdu,whdu,ra,dec,size,firsthdu=firsthdu,
        #             overlay_cat=ots,overlay_scale=scale,coords_color='red',
        #             coords_ra=r['RA'],coords_dec=r['DEC'],coords_lw=3,lw=2,
        #             save_name=wiseimage,no_labels=True,marker_ra=marker_ra,
        #             marker_dec=marker_dec,marker_lw=3,marker_color='cyan',
        #             title=title,peak=peak,lw_ellipse=3,ellipse_style=ls,
        #             ellipse_color='cyan',drlimit=1000,rms_use=rms,
        #             lofarlevel=2.5,noisethresh=0)

        
        with open(manifestname,'w') as manifest:
            manifest.write('%i,%s,%s,%s,%s,%f,%f,%f\n' % (i,psimage,csimage,pspimage,sourcename,ra,dec,size*3600.0))

        os.system('mogrify -quality 90 -trim '+sourcename+'*.png')

