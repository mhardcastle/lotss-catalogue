#!/usr/bin/python

'''
get_intersecting_sources
find which sources in a catalogue intersect, enclose or are enclosed by another source...
'''

from lofar_source_sorter import Mask, Masks_disjoint_complete
import numpy as np
from shapely.geometry.polygon import LinearRing
from shapely.geometry import Polygon
from astropy.table import Table
import astropy.coordinates as ac
import astropy.units as u
import os


path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'
lofarcat_file_srt = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.srl.sorted.fits'






lofarcat = Table.read(lofarcat_file_srt)


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
    mp = ea.intersection(eb)
    #encloses = ea.contains(eb)
    #encloses = ea.contains(Point(ellipse2[0], ellipse2[1]))
    pa = Polygon(a)
    pb = Polygon(b)
    encloses = pa.contains(pb)

    x = [p.x for p in mp]
    y = [p.y for p in mp]
    if len(x) > 0:
        intersects = True
    else:
        intersects =  False
    
    return intersects, encloses


calculate_intersections = True
if calculate_intersections:
    c = ac.SkyCoord(lofarcat['RA'], lofarcat['DEC'], unit="deg")
    cintersects = np.zeros(len(lofarcat), dtype=bool)
    cencloses = np.zeros(len(lofarcat), dtype=bool)
    cenclosed = np.zeros(len(lofarcat), dtype=bool)
    for i in range(len(lofarcat)):
        #print i
        l = lofarcat[i]
        if l['NN_sep'] > l['Maj']:  continue
        
        ci = ac.SkyCoord([l['RA']], [l['DEC']], unit='deg')
        idx, idxself, sep, _ = ci.search_around_sky(c, 2*l['Maj']*u.arcsec)
        if len(idx) == 1:   
            continue
        
        ii = idx
        
        for j in ii:
            if j == i: continue
            l1 = lofarcat[j]
            
            intersects, encloses = ellipses_intersect( (l['RA'], l['DEC'], l['Maj']/2., l['Min']/2, l['PA']), (l1['RA'], l1['DEC'], l1['Maj']/2., l1['Min']/2, l1['PA']) )
            
            print i, j, intersects, encloses, (l['RA'], l['DEC'], l['Maj']/2., l['Min']/2, l['PA']), (l1['RA'], l1['DEC'], l1['Maj']/2., l1['Min']/2, l1['PA'])
            
            if intersects:
                cintersects[i] = True
                cintersects[j] = True
            if encloses:
                cencloses[i] = True
                cenclosed[j] = True
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
    



## write output file

if os.path.exists(lofarcat_file_srt):
    os.remove(lofarcat_file_srt)
lofarcat.write(lofarcat_file_srt)