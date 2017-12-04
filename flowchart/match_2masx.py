#!/usr/bin/python

'''
match_2masx
find which sources potentially match to a 2MASX source
'''

from lofar_source_sorter import Mask, Masks_disjoint_complete
import numpy as np
from matplotlib import patches
from astropy.table import Table, join, Column
import astropy.coordinates as ac
import astropy.units as u
import os


path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'
lofarcat_file = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.srl.fixed.fits'
lofarcat_file_psrt = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.srl.fixed.presort.fits'

lofarcat = Table.read(lofarcat_file)


## get 2MASX information (from 'fixed' catalgoue)
xsc_file = path+'2MASX_hetdex_fix.fits'
xsc = Table.read(xsc_file)

Nxsc = len(xsc)

c = ac.SkyCoord(lofarcat['RA'], lofarcat['DEC'], unit="deg")
cxsc = ac.SkyCoord(xsc['ra'], xsc['dec'], unit="deg")
f_nn_idx,f_nn_sep2d,f_nn_dist3d = ac.match_coordinates_sky(c,cxsc,nthneighbor=1)

xsc_nn = xsc[f_nn_idx]




def accept_match_2mass(mask, lcat, xcat, plot=False, selname=None):

    
    if plot:
        f,ax = pp.paper_single_ax()
    if selname is not None:
        nn = lcat['Source_Name']==selname
        mask = mask & nn
        
    idx = np.arange(len(lcat))
    iinds = idx[mask]
    inellipse = np.zeros(len(lcat ), dtype=bool)
    for i in iinds:
        
        tl = lcat[i]
        tx = xcat[i]
        
        lerr = np.max((tl['E_RA'],tl['E_DEC']))
            
        # assumning flat sky here...
        g_ell_center = (tx['ra'], tx['dec'])
        r_a = (tx['r_ext']+ lerr )  / 3600. #to deg
        r_b = r_a*tx['k_ba']
        angle = 90.-tx['k_phi']  #(anticlockwise from x-axis)
        rangle = angle *np.pi/180.

        cos_angle = np.cos(np.radians(180.-angle))
        sin_angle = np.sin(np.radians(180.-angle))

        xc = tl['RA'] - g_ell_center[0]
        yc = tl['DEC'] - g_ell_center[1]

        xct = xc * cos_angle - yc * sin_angle
        yct = xc * sin_angle + yc * cos_angle 

        rad_cc = (xct/r_a)**2. + (yct/r_b)**2.
        
        if rad_cc <= 1:
            inellipse[i] = 1
            
        if plot:
            g_ellipse = patches.Ellipse(g_ell_center, 2*r_a, 2*r_b, angle=angle, fill=False, ec='green', linewidth=2)

            ax.plot(g_ell_center[0], g_ell_center[1], 'g.')
            ax.plot(tl['RA'], tl['DEC'], 'k.')
            ax.add_patch(g_ellipse)
            
            l_ellipse = patches.Ellipse((tl['RA'], tl['DEC']), 2*lerr/ 3600., 2*tl['Min']/ 3600., angle=90.-tl['PA'], fill=False, ec='blue', linewidth=2)
            ax.add_patch(l_ellipse)
            
            if rad_cc <= 1:
                ax.plot(tl['RA'], tl['DEC'], 'r+')
            
            mell = g_ellipse.contains_point((tl['RA'], tl['DEC']), radius=lerr/3600)
            if mell:
                ax.plot(tl['RA'], tl['DEC'], 'gx')
            
            ax.plot(g_ell_center[0], g_ell_center[1], 'g+')
            ax.plot(tl['RA'], tl['DEC'], 'b.')
    if plot:
        ax.set_xlabel('RA')
        ax.set_ylabel('DEC')
        ax.invert_xaxis()
        ax.axis('equal')
        
    return inellipse
        
#xmatch0 = f_nn_sep2d.value*u.deg < np.array(xsc_nn['r_ext'])*u.arcsec
xmatch1 = f_nn_sep2d.value*u.deg < np.array(xsc_nn['r_ext'] + np.max((lofarcat['E_RA'],lofarcat['E_DEC']),axis=0))*u.arcsec

inellipse = accept_match_2mass(xmatch1, lofarcat, xsc_nn)
xmatch = xmatch1 & inellipse

Xhuge =  xmatch & (xsc_nn['r_ext'] >= 60.)
Xsmall =  xmatch & (xsc_nn['r_ext'] >= 0.) & (xsc_nn['r_ext'] < 60.)

# add the columns if we've not yet run this script
if '2MASX' not in lofarcat.colnames:
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=bool),'2MASX'))
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype='S20'),'2MASX_name'))
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=float),'2MASX_ra'))
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=float),'2MASX_dec'))
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=float),'2MASX_size'))
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=bool),'2MASX_match_large'))
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=bool),'2MASX_match'))
    
for m in [Xhuge, Xsmall]:
    lofarcat['2MASX'][m]  = m[m]
    lofarcat['2MASX_name'][m]  = xsc_nn['designation'][m]
    lofarcat['2MASX_ra'][m]  = xsc_nn['ra'][m]
    lofarcat['2MASX_dec'][m]  = xsc_nn['dec'][m]
    lofarcat['2MASX_size'][m]  = xsc_nn['r_ext'][m]
    
    
lofarcat['2MASX_match_large'] = Xhuge
lofarcat['2MASX_match'] = Xsmall



print '{n:n} lofar sources have a possible 2MASX match'.format(n=sum(lofarcat['2MASX']))
print '{n:n} 2MASX sources out of {nall:n} have a possible LOFAR match'.format(n=len(np.unique(lofarcat['2MASX_name']))-1,nall=Nxsc)


## write output file

if os.path.exists(lofarcat_file_psrt):
    os.remove(lofarcat_file_psrt)
lofarcat.write(lofarcat_file_psrt)
