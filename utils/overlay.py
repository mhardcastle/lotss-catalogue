# overlay code ripped out of show_overlay's various incarnations

# the idea is to separate the plotting of the overlay from the
# determination of what maps to use

# lofarhdu must be a flattened (2D) LOFAR map. opthdu can be any type
# of optical map. rms_use must be the local rms

from __future__ import print_function
import os
from time import sleep
import numpy as np
import aplpy
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.table import Table
import random
import string

def choices(string, k=1):
    outstr=''
    for i in range(k):
        outstr+=string[np.random.randint(len(string))]
    return outstr

def tempnam():
    return ''.join(choices(string.ascii_uppercase + string.digits, k=10))

def find_noise(a):
    b=a.flatten()
    for i in range(10):
        m=np.nanmean(b)
        s=np.nanstd(b)
        b=b[b<(m+5.0*s)]
    return m,s

def find_noise_area(hdu,ra,dec,size,channel=0,true_max=False,debug=False):
    # ra, dec, size in degrees
    size/=1.5
    if len(hdu[0].data.shape)==2:
        cube=False
        ysize,xsize=hdu[0].data.shape
    else:
        channels,ysize,xsize=hdu[0].data.shape
        cube=True
    w=WCS(hdu[0].header)
    ras=[ra-size,ra-size,ra+size,ra+size]
    decs=[dec-size,dec+size,dec-size,dec+size]
    xv=[]
    yv=[]
    for r,d in zip(ras,decs):
        if cube:
            x,y,c=w.wcs_world2pix(r,d,0,0)
        else:
            x,y=w.wcs_world2pix(r,d,0)
        xv.append(x)
        yv.append(y)
    if debug: print(xv,yv)
    xmin=int(min(xv))
    if xmin<0: xmin=0
    xmax=int(max(xv))
    if xmax>=xsize: xmax=xsize-1
    ymin=int(min(yv))
    if ymin<0: ymin=0
    ymax=int(max(yv))
    if ymax>=ysize: ymax=ysize-1
    if debug: print(xmin,xmax,ymin,ymax)
    if cube:
        subim=hdu[0].data[channel,ymin:ymax,xmin:xmax]
    else:
        subim=hdu[0].data[ymin:ymax,xmin:xmax]
    if debug: subim.shape
    mean,noise=find_noise(subim)
    if true_max:
        vmax=np.nanmax(subim)
    else:
        for i in range(5,0,-1):
            vmax=np.nanmedian(subim[subim>(mean+i*noise)])
            if not(np.isnan(vmax)):
                break
    return mean,noise,vmax

def show_overlay(lofarhdu,opthdu,ra,dec,size,firsthdu=None,vlasshdu=None,rms_use=None,bmaj=None,bmin=None,bpa=None,title=None,save_name=None,plotpos=None,ppsize=750,block=True,interactive=False,plot_coords=True,overlay_cat=None,lw=1.0,show_lofar=True,no_labels=False,show_grid=True,overlay_region=None,overlay_scale=1.0,circle_radius=None,coords_color='white',coords_lw=1,coords_ra=None,coords_dec=None,marker_ra=None,marker_dec=None,marker_color='white',marker_lw=3,noisethresh=1,lofarlevel=2.0,contour_color='yellow',first_color='lightgreen',vlass_color='salmon',drlimit=500,interactive_handler=None,peak=None,minimum=None,ellipse_color='red',lw_ellipse=3,ellipse_style='solid',logfile=None,sourcename=None,vmax_cap=None,cmap='jet',lofar_colorscale=False,rotate_north=True,figure=None,figsize=None):
    '''
    show_overlay: make an overlay using AplPY.
    lofarhdu: the LOFAR cutout to use for contours
    opthdu: the optical image to use (may be 2D image or RGB cube)
    ra, dec: position to use, in degrees
    size: size in arcseconds
    firsthdu: FIRST cutout to use for contours or None if not required
    vlasshdu: VLASS cutout to use for contours or None if not required
    rms_use: Use a particular LOFAR rms or None to compute it from cutout
    bmaj, bmin, bpa: Beam sizes to use or None to take from the headers
    title: title at the top of the plot or None for no title
    save_name: filename to save image as
    plotpos: list of RA, Dec to plot as markers
    ppsize: size of the markers
    block: if True, wait before exiting if not saving an image
    interactive: if True, allow interaction with the image
    plot_coords: plot a set of co-ordinates as a large cross
    overlay_cat: plot a catalogue (astropy Table) as a set of ellipses with BMaj, BMin, BPA, or None to disable
    lw: line width for contours in Matplotlib units
    show_lofar: if True, show the LOFAR contours
    no_labels: if True, don't produce co-ordinate labels
    show_grid: if True, show a co-ordinate grid
    overlay_region: overlay a ds9-style region if not None
    overlay_scale: scaling factor for catalogue in overlay_cat
    circle_radius: if not None, draw a circle of this radius around the image centre
    coords_color: colour for plot_coords
    coords_lw: line width for plot_coords
    coords_ra, coords_dec: use these co-ordinates for plot_coords
    marker_ra, marker_dec: old marker list specification; deprecated, use plotpos instead
    marker_color,marker_lw: attributes for markers
    noisethresh: scaling noise for the optical images (lower cutoff at mean + rms*noisethresh)
    lofarlevel: lowest LOFAR contour (units sigma)
    first_color: colour for FIRST contours
    vlass_color: colour for VLASS contours
    drlimit: dynamic range limit for LOFAR data, default 500
    interactive_handler: function to call for clicks on an interactive map
    peak: use this LOFAR flux density as the peak value; if None, calculate from the map
    minimum: use this LOFAR flux density as the minimum contour level, over-riding DR limitations; if None, only use drlimit or lofarlevel
    ellipse_color: colour for the ellipses in overlay_cat: either a string or a list with entries for each ellipse
    lw_ellipse: line width for the ellipses; either a string or a list with entries for each ellipse
    ellipse_style: line style for the ellipses; either a string or a list with entries for each ellipse

    '''
    kwargs={}
    if figsize is not None:
        kwargs['figsize']=figsize
    if lofarhdu is None:
        print('LOFAR HDU is missing, not showing it')
        show_lofar=False
    try:
        from matplotlib.cbook import MatplotlibDeprecationWarning
        import warnings
        warnings.simplefilter('ignore', MatplotlibDeprecationWarning)
    except:
        print('Cannot hide warnings')

    print('===== Doing overlay at',ra,dec,'with size',size,'=====')
    print('===== Title is',title,'=====')
    
    if coords_ra is None:
        coords_ra=ra
        coords_dec=dec

    if show_lofar or lofar_colorscale:
        if peak is None:
            lofarmax=np.nanmax(lofarhdu[0].data)
            print('Local maximum flux density is',lofarmax)
        else:
            print('Using user-specified peak flux of',peak)
            lofarmax=peak
        if rms_use is None:
            rms_use=find_noise_area(lofarhdu,ra,dec,size)[1]
            print('Using LOFAR rms',rms_use)
        dr_minimum=lofarmax/drlimit
        print('Dynamic range minimum flux density is',dr_minimum,'with DR limit',drlimit)
        rms_minimum=rms_use*lofarlevel
        print('RMS range minimum flux density is',rms_minimum,'at',lofarlevel,'sigma')
        if minimum<dr_minimum:
            dr_minimum=minimum # override e.g. if there is a faint source near a bright one
            
        minlevel=max([dr_minimum,rms_minimum])
        print('Using minimum level',minlevel)
        if peak is not None and peak/minlevel<10:
            lofarmax=np.nanmax(lofarhdu[0].data) # revert
        levels=minlevel*2.0**np.linspace(0,14,30)

    rgbname=None
    hdu=opthdu

    if len(hdu[0].data.shape)>2:
        # This is an RGB image. Try to do something sensible with it
        vmins=[]
        vmaxes=[]
        for i in range(3):
            mean,noise,vmax=find_noise_area(hdu,ra,dec,size,channel=i)
            vmin=mean+noisethresh*noise
            print('channel',i,mean,noise,vmin,vmax)
            vmins.append(vmin)
            vmaxes.append(vmax)
        vmax=np.nanmax(vmaxes)
        if vmax_cap is not None and vmax>vmax_cap:
            vmax=vmax_cap
        rgbname='rgb_'+tempnam()+'.png'
        aplpy.make_rgb_image(hdu.filename(),rgbname,stretch_r='log',stretch_g='log',stretch_b='log',vmin_r=vmins[0],vmin_g=vmins[1],vmin_b=vmins[2],vmax_r=vmax,vmax_g=vmax,vmax_b=vmax)
        f=aplpy.FITSFigure(rgbname,north=rotate_north,**kwargs)
        print('centring on',ra,dec,size)
        f.recenter(ra,dec,width=size,height=size)
        f.show_rgb()
        if logfile is not None:
            mytuple=tuple([sourcename]+vmins+vmaxes+[vmax])
            print(mytuple)
            logfile.write('%s %f %f %f %f %f %f %f\n' % mytuple)
    elif lofar_colorscale:
        print('Doing a LOFAR colour scale')
        f = aplpy.FITSFigure(hdu,north=rotate_north,**kwargs)
        print('centring on',ra,dec,size)
        f.recenter(ra,dec,width=size,height=size)
        f.show_colorscale(vmin=rms_use,vmax=lofarmax,stretch='log',cmap=cmap)
   
    else:
        f = aplpy.FITSFigure(hdu,north=rotate_north,**kwargs)
        print('centring on',ra,dec,size)
        f.recenter(ra,dec,width=size,height=size)
        mean,noise,vmax=find_noise_area(hdu,ra,dec,size)
        if not np.isnan(mean) and not np.isnan(vmax) and not np.isnan(noise):
            print('Optical parameters are',mean,noise,vmax)
            f.show_colorscale(vmin=mean+noisethresh*noise,vmax=vmax,stretch='log',cmap=cmap)
        else:
            print('*** Warning -- cannot show optical image, nans present ***')
        #f.show_colorscale(stretch='log')

    if bmaj is not None:
        f._header['BMAJ']=bmaj
        f._header['BMIN']=bmin
        f._header['BPA']=bpa

    if show_lofar: f.show_contour(lofarhdu,colors=contour_color,linewidths=lw, levels=levels)

    if firsthdu is not None:
        firstrms=find_noise_area(firsthdu,ra,dec,size)[1]
        print('Using FIRST rms',firstrms)
        firstlevels=firstrms*3*2.0**np.linspace(0,14,30)
        f.show_contour(firsthdu,colors=first_color,linewidths=lw, levels=firstlevels)

    if vlasshdu is not None:
        vlassrms=find_noise_area(vlasshdu,ra,dec,size)[1]
        print('Using VLASS rms',vlassrms)
        vlasslevels=vlassrms*3*2.0**np.linspace(0,14,30)
        f.show_contour(vlasshdu,colors=vlass_color,linewidths=lw, levels=vlasslevels)

    if bmaj is not None:
        f.add_beam()
        f.beam.set_corner('bottom left')
        f.beam.set_edgecolor('red')
        f.beam.set_facecolor('white')

    if plot_coords:
        f.show_markers(coords_ra,coords_dec,marker='+',facecolor=coords_color,edgecolor=coords_color,linewidth=coords_lw,s=1500,zorder=100)

    if marker_ra is not None:
        f.show_markers(marker_ra,marker_dec,marker='x',facecolor=marker_color,edgecolor=marker_color,linewidth=marker_lw,s=1500,zorder=140)
        

    if plotpos is not None:
        if not isinstance(plotpos,list):
            plotpos=[(plotpos,'x'),]
        for element in plotpos:
            if len(element)==4:
                # 4 elements: specify edgecolor and fill color separately, can be none
                t,marker,edgecolor,facecolor=element
            elif len(element)==3:
                # 3 elements, face color is edge color
                t,marker,edgecolor=element
                facecolor=edgecolor # was 'None'
            else:
                # 2 elements, mark in white
                t,marker=element
                edgecolor='white'
                facecolor=edgecolor
            if t is not None and len(t)>0:
                f.show_markers(t['ra'],t['dec'],marker=marker,facecolor=facecolor,edgecolor=edgecolor,linewidth=2,s=ppsize,zorder=100)

    if circle_radius is not None:
        f.show_circles([ra,],[dec,],[circle_radius,],facecolor='none',edgecolor='cyan',linewidth=5)

    if overlay_cat is not None:
        t=overlay_cat
        f.show_ellipses(t['RA'],t['DEC'],t['Maj']*2/overlay_scale,t['Min']*2/overlay_scale,angle=90+t['PA'],edgecolor=ellipse_color,linewidth=lw_ellipse,linestyle=ellipse_style,zorder=100)

    if overlay_region is not None:
        f.show_regions(overlay_region)

    if no_labels:
        f.axis_labels.hide()
        f.tick_labels.hide()
    if show_grid:
        f.add_grid()
        f.grid.show()
        f.grid.set_xspacing(1.0/60.0)
        f.grid.set_yspacing(1.0/60.0)

    if title is not None:
        plt.title(title)
    plt.tight_layout()
    if save_name is None:
        plt.show(block=block)
    else:
        plt.savefig(save_name)
        plt.close()

    def onclick(event):
        xp=event.xdata
        yp=event.ydata
        xw,yw=f.pixel2world(xp,yp)
        if event.button==2:
            print(title,xw,yw)

    if interactive:
        fig=plt.gcf()
        if interactive_handler is None:
            cid = fig.canvas.mpl_connect('button_press_event', onclick)
        else:
            I=interactive_handler(f)
            fig.canvas.mpl_connect('button_press_event', I.onclick)

    if rgbname is not None:
        os.remove(rgbname)
        
    return f
