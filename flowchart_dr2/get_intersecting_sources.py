#!/usr/bin/python

'''
get_intersecting_sources
find which sources in a catalogue intersect, enclose or are enclosed by another source...
'''

from lofar_source_sorter_dr2 import Mask, Masks_disjoint_complete
import numpy as np
from shapely.geometry.polygon import LinearRing
from shapely.geometry import Polygon
from astropy.table import Table, Column
import astropy.coordinates as ac
import astropy.units as u
import os

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import utils.plotting as pp

#path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'
#lofarcat_file_srt = path+'LOFAR_HBA_T1_DR1_catalog_v0.95_masked.srl.fixed.sorted.fits'




path = '/data2/wwilliams/projects/lofar_surveys/LoTSS-DR2-Feb2020/'
lofarcat_file_srt = path+'LoTSS_DR2_rolling.srl_0h.lr.presort.fits'




lofarcat = Table.read(lofarcat_file_srt)

# temporary selection for testing
#ti = (lofarcat['RA']>330) & (lofarcat['RA']<332) & (lofarcat['DEC']>25) & (lofarcat['DEC']<27) 
#ti = (lofarcat['RA']>331.52) & (lofarcat['RA']<331.7) & (lofarcat['DEC']>25.27) & (lofarcat['DEC']<25.39) 
#lofarcat = lofarcat[ti]

def ellipse_polyline(ellipses, n=100):
    t = np.linspace(0, 2*np.pi, n, endpoint=True)
    st = np.sin(t)
    ct = np.cos(t)
    result = []
    for x0, y0, a, b, angle in ellipses:
        angle = np.deg2rad(angle)
        sa = np.sin(angle)
        ca = np.cos(angle)
        p = np.empty((n, 2))
        p[:, 0] = x0 + a * ca * ct - b * sa * st
        p[:, 1] = y0 + a * sa * ct + b * ca * st
        result.append(p)
    return result

def ellipses_intersect(ellipse1, ellipse2):
    a, b = ellipse_polyline((ellipse1, ellipse2))
    ea = LinearRing(a)
    eb = LinearRing(b)
    #mp = ea.intersection(eb)
    ##encloses = ea.contains(eb)
    ##encloses = ea.contains(Point(ellipse2[0], ellipse2[1]))
    pa = Polygon(a)
    pb = Polygon(b)
    #encloses = pa.contains(pb)

    #c = mp.coords
    ##x = [p.x for p in ]
    ##y = [p.y for p in mp.coords]
    #if len(c) > 0:
        #intersects = True
    #else:
        #intersects =  False
    
    return ea.crosses(eb), pa.contains(pb), pb.contains(pa)

# work in degrees
lofarcat['Maj'] = lofarcat['Maj'] /3600.
lofarcat['Min'] = lofarcat['Min'] /3600.
lofarcat['NN_Maj'] = lofarcat['NN_Maj'] /3600.
lofarcat['NN_sep'] = lofarcat['NN_sep'] /3600.



calculate_intersections = True
if calculate_intersections:
    c = ac.SkyCoord(lofarcat['RA'], lofarcat['DEC'], unit="deg")
    cintersects = np.zeros(len(lofarcat), dtype=bool)
    cencloses = np.zeros(len(lofarcat), dtype=bool)
    cenclosed = np.zeros(len(lofarcat), dtype=bool)
    
    ells  = []
    ellcols = []
    for i in range(len(lofarcat)):
        #print i
        l = lofarcat[i]
        
        #if l['NN_sep'] > l['Maj']/10.:  continue
        
        ## this may miss sources where the NN is not the overlapper? is this possible?
        # potential overlap sources lie within a separation of the sizes of the two sources
        if l['NN_sep'] > (l['Maj']+l['NN_Maj']):  continue
        
        ci = ac.SkyCoord([l['RA']], [l['DEC']], unit='deg')
        idx, idxself, sep, _ = ci.search_around_sky(c, (l['Maj']+l['NN_Maj'])*u.deg)
        if len(idx) == 1:   
            continue
        
        ii = idx
        
        # all in gray
        col = 'gray'
        
        
        #f1,ax1 = pp.paper_single_ax()
        #ax1.set_title(l['Source_Name'])
        
        for j in ii:
            if j == i: continue
            l1 = lofarcat[j]
            
            intersects, encloses, enclosed = ellipses_intersect( (l['RA'], l['DEC'], l['Maj']/2., l['Min']/2, l['PA']+90), (l1['RA'], l1['DEC'], l1['Maj']/2., l1['Min']/2, l1['PA']+90) )
            ellipse1, ellipse2 = (l['RA'], l['DEC'], l['Maj']/2., l['Min']/2, l['PA']+90), (l1['RA'], l1['DEC'], l1['Maj']/2., l1['Min']/2, l1['PA']+90)
            ei, ej = ellipse_polyline((ellipse1, ellipse2) )
        
            
            print(i, j, len(ii), intersects, encloses, enclosed, (l['RA'], l['DEC'], l['Maj']/2., l['Min']/2, l['PA']), (l1['RA'], l1['DEC'], l1['Maj']/2., l1['Min']/2, l1['PA']))
            
            if intersects:
                cintersects[i] = True
                #cintersects[j] = True
                col = 'C1'
            if encloses:
                cencloses[i] = True
                col = 'C0'
            if enclosed:
                cenclosed[i] = True
                col = 'C2'
                
            #ax1.plot(ei.T[0], ei.T[1], c=col)    
            #ax1.plot(ej.T[0], ej.T[1], c='gray')

        #plt.show(f1)
        #cc = input('cont')
        #plt.close(f1)
            
        #ells.append(Ellipse(xy=(l['RA'], l['DEC']), width=l['Maj']/3600., height=l['Min']/3600., angle=l['PA']))
        ells.append((l['RA'], l['DEC'], l['Maj'], l['Min'], l['PA']))
        ellcols.append(col)
        
    if 'Encloses' not in lofarcat.colnames:
        lofarcat.add_column(Column(cencloses, 'Encloses'))
    else:
        lofarcat['Encloses'] = cencloses
    if 'Enclosed' not in lofarcat.colnames:
        lofarcat.add_column(Column(cenclosed, 'Enclosed'))
    else:
        lofarcat['Enclosed'] = cenclosed
    if 'Intersects' not in lofarcat.colnames:
        lofarcat.add_column(Column(cintersects, 'Intersects'))
    else:
        lofarcat['Intersects'] = cintersects
        

    M_encloses = Mask(lofarcat['Encloses'],
                    'Encloses',
                    'encloses',
                    qlabel='encloses?')
    M_enclosed = Mask(lofarcat['Enclosed'],
                    'Enclosed',
                    'enclosed',
                    qlabel='enclosed?')
    M_intersects = Mask(lofarcat['Intersects'],
                    'Intersects',
                    'intersects',
                    qlabel='intersects?')
    M_encloses.make_sample(lofarcat)
    M_enclosed.make_sample(lofarcat)
    M_intersects.make_sample(lofarcat)
    

plot = False
if plot:
    f,ax = pp.paper_single_ax(TW=12)

    for ell, ellcol in zip(ells,ellcols):
        li = Ellipse(xy=(ell[0], ell[1]), width=ell[2], height=ell[3], angle=ell[4]+90.)
        li.set_edgecolor(ellcol)
        #l.set_facecolor(None)
        li.set_facecolor('none')
        ax.add_artist(li)
        
    ax.scatter(lofarcat['RA'], lofarcat['DEC'],s=1,c='gray')
    ax.set_ylim(lofarcat['DEC'].min(), lofarcat['DEC'].max())
    ax.set_xlim(lofarcat['RA'].max(), lofarcat['RA'].min())


## write output file

if os.path.exists(lofarcat_file_srt):
    os.remove(lofarcat_file_srt)
lofarcat.write(lofarcat_file_srt)
