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
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import utils.plotting as pp



#from multiprocessing import Process, Value, Array, shared_memory, Pool
from multiprocessing import Pool, Array, Value

#def f(n, a):
    #n.value = 3.1415927
    #for i in range(len(a)):
        #a[i] = -a[i]

#def test_multi():
    #num = Value('d', 0.0)
    #arr = Array('i', range(10))
    #arr2 = Array('i', range(10))

    #p = Process(target=f, args=(num, arr))
    #p.start()
    #p.join()

    
    

#path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'
#lofarcat_file_srt = path+'LOFAR_HBA_T1_DR1_catalog_v0.95_masked.srl.fixed.sorted.fits'


if len(sys.argv) == 1:
    print("Usage is : python get_intersecting_sources.py field_code ")
    print('E.g.: python get_intersecting_sources.py 0 ')
    sys.exit(1)

h = str(sys.argv[1])
if 'h' not in h:
    h+='h'
if h not in  ['0h','13h','n0h','n13h','s0h','s13h']:
    print('unknown field code (should be 0h or 13h)',h)
    if 'subcat' not in h:
        print('not a subcat')
        sys.exit(1)

path = '/Users/w.williams/projects/lofar_surveys/DR2/'
lofarcat_file_srt = path+'LoTSS_DR2_{version}.srl_{h}.lr-full.presort.fits'.format(version=version,h=h)




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




calculate_intersections = True
if calculate_intersections:
    #c = ac.SkyCoord(lofarcat['RA'], lofarcat['DEC'], unit="deg")
    #cintersects = np.zeros(len(lofarcat), dtype=bool)
    #cencloses = np.zeros(len(lofarcat), dtype=bool)
    #cenclosed = np.zeros(len(lofarcat), dtype=bool)
    
    

    def get_intersections(i):
        
        isq = (abs(np.array(Ara[:]) -Ara[i]) < 2*MaxSize.value) & (abs(np.array(Adec[:]) -Adec[i]) < 2*MaxSize.value)
        #print(i, np.sum(isq))
        
        if np.sum(isq) == 0:
            print(i,'isol')
            return False, False, False
        #print(i,np.sum(isq))
        
          
        ci = ac.SkyCoord([Ara[i]], [Adec[i]], unit='deg')
        
        isq = np.where(isq)[0]
        
        ## this may miss sources where the NN is not the overlapper? is this possible?
        # potential overlap sources lie within a separation of the sizes of the two sources
        #if Annsep[i] > (Amaj[i]+Annmaj[i]):  
            #return False, False, False
        
        cintersects = False
        cencloses = False
        cenclosed = False
        
        for j in isq:
            if j == i: continue
        
            c = ac.SkyCoord(Ara[j], Adec[j], unit="deg")  
            
            # rough check on shapes
            if  ci.separation(c) >= (Amaj[i]+Amaj[j])*u.deg:
                continue
            
            
            intersects, encloses, enclosed = ellipses_intersect( (Ara[i], Adec[i], Amaj[i]/2., Amin[i]/2, Apa[i]+90), (Ara[j], Adec[j], Amaj[j]/2., Amin[j]/2, Apa[j]+90) )
            #ellipse1, ellipse2 = (Ara[i], Adec[i], Amaj[i]/2., Amin[i]/2, Apa[i]+90), (Ara[j], Adec[j], Amaj[j]/2., Amin[j]/2, Apa[j]+90)
            #ei, ej = ellipse_polyline((ellipse1, ellipse2) )
        
            
            #print(i, j, len(idx), intersects, encloses, enclosed, (Ara[i], Adec[i], Amaj[i]/2., Amin[i]/2, Apa[i]), (Ara[j], Adec[j], Amaj[j]/2., Amin[j]/2, Apa[j]))
            
            if intersects:
                cintersects = True
            if encloses:
                cencloses = True
            if enclosed:
                cenclosed = True
            
        #print(i,cintersects,cencloses,cenclosed)
                
        return cintersects, cencloses, cenclosed

    #lofarcat = Table.read('/home/wendy/projects/lba_bootes_deep/bootes_deep_lba.cat.fits')
    #lofarcat = Table.read('/Users/w.williams/projects/lofar_surveys/DR2/LoTSS_DR2_{version}.srl_13h.lr-full.presort.fits')


    # work in degrees
    lofarcat['Maj'] = lofarcat['Maj'] /3600.
    lofarcat['Min'] = lofarcat['Min'] /3600.
    lofarcat['NN_Maj'] = lofarcat['NN_Maj'] /3600.
    lofarcat['NN_sep'] = lofarcat['NN_sep'] /3600.

    # maked shared arrays - saves memory on multiproc
    #MaxSize = Value('d', lofarcat['Maj'].max(), lock=False)
    MaxMaj = lofarcat['Maj'].max()
    MaxMajSmall = 60./3600
    
    
    ind_large = np.where(lofarcat['Maj'] > MaxMajSmall)[0]
    ind_small = np.where(lofarcat['Maj'] <= MaxMajSmall)[0]
    # get all the indicies for sources near large sources
    ind_near_large = np.array([],dtype=int)
    for i in ind_large:
        #print(i)
        new_ind = np.where((np.abs(lofarcat['RA'][i] - lofarcat['RA']) < 2*MaxMaj) & (np.abs(lofarcat['DEC'][i] - lofarcat['DEC']) < 2*MaxMaj))[0]
        ind_near_large = np.hstack((ind_near_large, new_ind))
    ind_near_large = np.unique(ind_near_large)
    
    MaxSize = Value('d', MaxMajSmall, lock=False) # first restrict to 60"
    
    Ara = Array('d', lofarcat['RA'], lock=False)
    Adec = Array('d', lofarcat['DEC'], lock=False)
    Amaj = Array('d', lofarcat['Maj'], lock=False)
    Amin = Array('d', lofarcat['Min'], lock=False)
    Apa = Array('d', lofarcat['PA'], lock=False)
    Annsep = Array('d', lofarcat['NN_sep'], lock=False)
    Annmaj = Array('d', lofarcat['NN_Maj'], lock=False)
    Ncat = len(lofarcat)
    Cintersects = np.zeros(Ncat,dtype=bool)
    Cencloses = np.zeros(Ncat,dtype=bool)
    Cenclosed = np.zeros(Ncat,dtype=bool)

    # should clear memory
    lofarcat = 1
    del lofarcat

    print('getting intersections', len(ind_small))
    with Pool(processes=12) as pool:
        result = pool.map(get_intersections, ind_small)
        #result = pool.map(get_intersections, range(100))
    result = np.array(result)
    Cintersects[ind_small] = result[:,0]
    Cencloses[ind_small] = result[:,1]
    Cenclosed[ind_small] = result[:,2]
    
    print(np.sum(Cintersects),' intersects')
    print(np.sum(Cencloses),' enclosing')
    print(np.sum(Cenclosed),' enclosed')
    
    MaxSize = Value('d', MaxMaj, lock=False)

    print('getting intersections - big sources', len(ind_near_large))
    with Pool(processes=12) as pool:
        result = pool.map(get_intersections, ind_near_large)
    result = np.array(result)
    Cintersects[ind_near_large] = result[:,0]
    Cencloses[ind_near_large] = result[:,1]
    Cenclosed[ind_near_large] = result[:,2]

    print(np.sum(Cintersects),' intersects')
    print(np.sum(Cencloses),' enclosing')
    print(np.sum(Cenclosed),' enclosed')
        
        
    lofarcat = Table.read(lofarcat_file_srt)
        
    if 'Encloses' not in lofarcat.colnames:
        lofarcat.add_column(Column(Cencloses, 'Encloses'))
    else:
        lofarcat['Encloses'] = Cencloses
    if 'Enclosed' not in lofarcat.colnames:
        lofarcat.add_column(Column(Cenclosed, 'Enclosed'))
    else:
        lofarcat['Enclosed'] = Cenclosed
    if 'Intersects' not in lofarcat.colnames:
        lofarcat.add_column(Column(Cintersects, 'Intersects'))
    else:
        lofarcat['Intersects'] = Cintersects
        

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


# put back in arcsec  - don't do this since we reopened the file
#lofarcat['Maj'] = lofarcat['Maj'] *3600.
#lofarcat['Min'] = lofarcat['Min'] *3600.
#lofarcat['NN_Maj'] = lofarcat['NN_Maj'] *3600.
#lofarcat['NN_sep'] = lofarcat['NN_sep'] *3600.


lofarcat.keep_columns(['Source_Name','block','Encloses','Enclosed','Intersects'])
## write output file

if os.path.exists(lofarcat_file_srt):
    os.remove(lofarcat_file_srt)
lofarcat.write(lofarcat_file_srt)
