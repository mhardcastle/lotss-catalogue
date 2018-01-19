#!/usr/bin/python

'''
get_candidate_doubles
find which sources are possible doubles
'''

from lofar_source_sorter import Mask, Masks_disjoint_complete
import numpy as np
from shapely.geometry.polygon import LinearRing
from shapely.geometry import Polygon
from shapely.ops import cascaded_union
from descartes import PolygonPatch
#from skimage.morphology import convex_hull_image
#from skimage.feature import peak_local_max
#from skimage.measure import label, regionprops
from astropy.table import Table, Column
import astropy.coordinates as ac
import astropy.units as u
import os

import matplotlib.pyplot as plt

import plot_util as pp

def ellipse(x0,y0,a,b,pa,n=200):
    a=a/2.
    b=b/2.
    theta=np.linspace(0,2*np.pi,n,endpoint=False)
    st = np.sin(theta)
    ct = np.cos(theta)
    th = np.deg2rad(90-pa)
    sa = np.sin(th)
    ca = np.cos(th)
    p = np.empty((n, 2))
    p[:, 0] = x0 + a * ca * ct - b * sa * st
    p[:, 1] = y0 + a * sa * ct + b * ca * st
    return Polygon(p)



lLR_thresh = 0.36

path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'
lofarcat_file_srt = path+'LOFAR_HBA_T1_DR1_catalog_v0.95_masked.srl.fixed.sorted.fits'






lofarcat = Table.read(lofarcat_file_srt)


lofarcat_small = lofarcat[lofarcat['Maj'] < 15.]
flux = lofarcat_small['Total_flux']

c = ac.SkyCoord(lofarcat_small['RA'], lofarcat_small['DEC'], unit="deg")

idx, idxpair, sep, _ = c.search_around_sky(c, seplimit=60*u.arcsec)
#remove matches to self
ex0 = (sep.value>1e-6)
idx = idx[ex0]
idxpair = idxpair[ex0]
sep = sep[ex0]

print '{:n} pairs'.format(len(idx))

#remove duplicates
dups = (flux[idx]/flux[idxpair] < 1.)
idx = idx[~dups]
idxpair = idxpair[~dups]
sep = sep[~dups]

print '{:n} unique pairs'.format(len(idx))
 
fratr = (flux[idx]/flux[idxpair]) < 10
frat = (flux[idx]/flux[idxpair]) < 4
dist = ((flux[idx]+flux[idxpair]) >= (50*(sep.to('arcsec').value/100)**2.))
distr = ((flux[idx]+flux[idxpair]) >= (45*(sep.to('arcsec').value/100)**2.))


ra = (flux[idx]*lofarcat_small['RA'][idx] + flux[idxpair]*lofarcat_small['RA'][idxpair]) / (flux[idx]+flux[idxpair])
dec = (flux[idx]*lofarcat_small['DEC'][idx] + flux[idxpair]*lofarcat_small['DEC'][idxpair]) / (flux[idx]+flux[idxpair])

#sep = np.sqrt((lofarcat_small['RA'][idx] -lofarcat_small['RA'][idxpair]) **2. + (lofarcat_small['DEC'][idx] -lofarcat_small['DEC'][idxpair]) **2.)

ra2 = (lofarcat_small['RA'][idx] + lofarcat_small['RA'][idxpair]) / (2.)
dec2 = (lofarcat_small['DEC'][idx] + lofarcat_small['DEC'][idxpair]) / (2.)


pa = c[idx].position_angle(c[idxpair])# +90.*u.deg
theta = 90*u.deg - pa

#dra = (lofarcat_small['RA'][idx] -lofarcat_small['RA'][idxpair])
#dra.unit=''
#ddec = (lofarcat_small['DEC'][idx] -lofarcat_small['DEC'][idxpair])
#ddec.unit=''
#theta = u.rad*np.arctan2(ddec,dra)
#pa = 90*u.deg - theta


# projected major axis along line connecting two sources
#pmaj1 = np.abs(lofarcat_small['Maj'][idx] * np.cos(pa - lofarcat_small['PA'][idx]))
#pmaj2 = np.abs(lofarcat_small['Maj'][idxpair] * np.cos(pa - lofarcat_small['PA'][idxpair]))
#pmin1 = np.abs(lofarcat_small['Min'][idx] * np.cos(pa - lofarcat_small['PA'][idx]))
#pmin2 = np.abs(lofarcat_small['Min'][idxpair] * np.cos(pa - lofarcat_small['PA'][idxpair]))
pmaj1 = np.abs(lofarcat_small['Maj'][idx] * np.cos(pa - lofarcat_small['PA'][idx].to('degree')))
pmaj2 = np.abs(lofarcat_small['Maj'][idxpair] * np.cos(pa - lofarcat_small['PA'][idxpair].to('degree')))
pmin1 = np.abs(lofarcat_small['Min'][idx] * np.cos(pa - lofarcat_small['PA'][idx].to('degree')))
pmin2 = np.abs(lofarcat_small['Min'][idxpair] * np.cos(pa - lofarcat_small['PA'][idxpair].to('degree')))


pmin = 0.2* np.max((pmaj1, pmaj2, pmin1, pmin2), axis=0)

pmaj = 0.2*sep

for i in range(10):
    f,ax = pp.paper_single_ax()
    ellist = []
    for ij in [idx[i], idxpair[i]]:
        print lofarcat_small['Source_Name'][ij]
        ax.plot(lofarcat_small['RA'][ij], lofarcat_small['DEC'][ij],'x')
        newp = ellipse(lofarcat_small['RA'][ij], lofarcat_small['DEC'][ij], lofarcat_small['Maj'][ij]/3600., lofarcat_small['Min'][ij]/3600., lofarcat_small['PA'][ij])
        print lofarcat_small['RA'][ij], lofarcat_small['DEC'][ij], lofarcat_small['Maj'][ij], lofarcat_small['Min'][ij], lofarcat_small['PA'][ij]
        patch = PolygonPatch(newp,alpha=0.5)
        ax.add_patch(patch)
        ellist.append(newp)
        
    ax.plot(ra2[i], dec2[i], 'x')
    #newp = ellipse(ra2[i], dec2[i], pmaj[i], pmin[i]/3600., pa[i].value)
    newp = ellipse(ra2[i], dec2[i], pmaj[i].degree, pmin[i]/3600., pa[i].degree)
    patch = PolygonPatch(newp,alpha=0.5)
    ax.add_patch(patch)
    ax.plot(ra[i], dec[i], 'x')
    ax.invert_xaxis()
    ax.axis('equal')
        


f,ax = pp.paper_single_ax()
ax.hist(lofarcat_small['E_Maj'], bins=100, histtype='step')
ax.hist(pmaj.to('arcsec'), bins=100, histtype='step')

f,ax = pp.paper_single_ax()
ax.hist(lofarcat_small['E_Min'], bins=100, histtype='step')
ax.hist(pmin, bins=100, histtype='step')


lrcrit = np.log10(1+lofarcat_small['LR']) > lLR_thresh

outtable = Table([Column(ra, 'RA'),
                  Column(dec, 'DEC'),
                  Column(pmaj.to('arcsec'), 'E_Maj'),
                  Column(pmin*u.arcsec, 'E_Min'),
                  Column(pa.to('deg'), 'PA'),
                  Column(lofarcat_small['Source_Name'][idx], 'Source_Name_1'),
                  Column(lofarcat_small['RA'][idx], 'RA_1'),
                  Column(lofarcat_small['DEC'][idx], 'DEC_1'),
                  Column(lofarcat_small['Maj'][idx], 'Maj_1'),
                  Column(lofarcat_small['Min'][idx], 'Min_1'),
                  Column(lofarcat_small['PA'][idx], 'PA_1'),
                  Column(flux[idx], 'Total_flux_1'),
                  Column(lofarcat_small['S_Code'][idx], 'S_Code_1'),
                  Column(lofarcat_small['Source_Name'][idxpair], 'Source_Name_2'),
                  Column(lofarcat_small['RA'][idxpair], 'RA_2'),
                  Column(lofarcat_small['DEC'][idxpair], 'DEC_2'),
                  Column(lofarcat_small['Maj'][idxpair], 'Maj_2'),
                  Column(lofarcat_small['Min'][idxpair], 'Min_2'),
                  Column(lofarcat_small['PA'][idxpair], 'PA_2'),
                  Column(flux[idxpair], 'Total_flux_2'),
                  Column(lofarcat_small['S_Code'][idxpair], 'S_Code_2'),
                  Column(sep, 'Separation'),
                  Column(np.sum([lofarcat_small['S_Code'][idxpair]=='S', lofarcat_small['S_Code'][idx]=='S'],axis=0), 'NS'),
                  Column(np.sum([lrcrit[idxpair], lrcrit[idx]],axis=0), 'NLR'),
                  Column(dist, 'dist'),
                  Column(distr, 'dist_r'),
                  Column(frat, 'Frat'),
                  Column(fratr, 'Frat_r'),])
outtable.add_column(Column(outtable['NS']==0, 's0'))
outtable.add_column(Column(outtable['NS']==1, 's1'))
outtable.add_column(Column(outtable['NS']==2, 's2'))
outtable.add_column(Column(outtable['NLR']==0, 'lr0'))
outtable.add_column(Column(outtable['NLR']==1, 'lr1'))
outtable.add_column(Column(outtable['NLR']==2, 'lr2'))

outtable.write('LOFAR_HBA_T1_DR1_catalog_v0.95_pairs.fits', overwrite=True)


### output from pepe
lrpaircat = Table.read('/local/wwilliams/projects/radio_imaging/lofar_surveys/source_class/doubles/lofar_pairs_pw.fits')

outtable = lrpaircat


#sys.exit()


plt.figure()
plt.scatter(ra,dec, s=4)
#for i, ii in zip(idx,idxpair):
    #plt.plot([lofarcat_small['RA'][i],lofarcat_small['RA'][ii]],[lofarcat_small['DEC'][i],lofarcat_small['DEC'][ii]], c='C0')

## both are  an S source
#ssource =  (lofarcat_small['S_Code'][idx] =='S') & ( lofarcat_small['S_Code'][idxpair] =='S')
## one or the other is an M source
#msource =  np.logical_xor((lofarcat_small['S_Code'][idx] !='S'), ( lofarcat_small['S_Code'][idxpair] !='S'))
## both are an M source
#msource2 =  (lofarcat_small['S_Code'][idx] !='S') & ( lofarcat_small['S_Code'][idxpair] !='S')

ssource = outtable['s2']
msource = outtable['s1']
msource2 = outtable['s0']


## neither has an lr match
#nlr =  (np.log10(1+lofarcat_small['LR'][idx]) <= lLR_thresh) & (np.log10(1+lofarcat_small['LR'][idxpair]) <= lLR_thresh)
## one or the other has an lr match
#lr = np.logical_xor((np.log10(1+lofarcat_small['LR'][idx]) > lLR_thresh),(np.log10(1+lofarcat_small['LR'][idxpair]) > lLR_thresh))
## both have an lr match
#lr2 =  (np.log10(1+lofarcat_small['LR'][idx]) > lLR_thresh) & (np.log10(1+lofarcat_small['LR'][idxpair]) > lLR_thresh)

nlr = outtable['lr0']
lr = outtable['lr1']
lr2 = outtable['lr2']





plt.figure()
plt.hist(np.log10(outtable['lr'][ssource]), range=(0,3), bins=100, histtype='step')
plt.hist(np.log10(outtable['lr'][ssource&nlr]), range=(0,3), bins=100, histtype='step', label='0')
plt.hist(np.log10(outtable['lr'][ssource&lr]), range=(0,3), bins=100, histtype='step', label='1')
plt.hist(np.log10(outtable['lr'][ssource&lr2]), range=(0,3), bins=100, histtype='step', label='2')
plt.legend()



print '{:n} unique pairs with both 1 S source'.format(np.sum(ssource))
print '{:n} unique pairs with only 1 M source'.format(np.sum(msource))
print '{:n} unique pairs with both M source'.format(np.sum(msource2))
print
print '{:n} unique pairs with neither 1 lr match'.format(np.sum(nlr))
print '{:n} unique pairs with only 1 lr match'.format(np.sum(lr))
print '{:n} unique pairs with both lr matches'.format(np.sum(lr2))
print

def printline(sstr,mask):
    n = np.sum(mask)
    S =  '{:n} unique pairs with {:s}'.format(n, sstr)
    print S
    
    outtable[mask].write('doubles/'+sstr.replace(' ','_')+'.fits', overwrite=True)
    
    
strs = np.array(['both M', 'one M', 'both S', 'both S and both lr', 'both S and one lr', 'both S and neither lr'])
masks = np.array([ msource2, msource, ssource, ssource*lr2, ssource*lr, ssource*nlr])
for sstr, mask in zip(strs, masks):
    printline(sstr, mask)
    printline(sstr+' and Fratio', mask & frat)
    printline(sstr+' and not Fratio', mask & (~frat))
    printline(sstr+' and Fratio relax', mask & fratr)
    printline(sstr+' and dist', mask & dist)
    printline(sstr+' and not dist', mask & (~dist))
    printline(sstr+' and dist relax', mask & distr)
    printline(sstr+' and Fratio and dist', mask & frat & dist)
    printline(sstr+' and not Fratio and not dist', mask & ~(frat & dist))
    printline(sstr+' and Fratio and dist relax', mask & fratr & distr)
    print 

plt.figure()
plt.semilogy()
plt.scatter(sep.to('arcsec'), flux[idx]+flux[idxpair], s=2, alpha=0.2)
plt.scatter(sep.to('arcsec')[frat], (flux[idx]+flux[idxpair])[frat], s=2, alpha=0.2)
plt.scatter(sep.to('arcsec')[dist], (flux[idx]+flux[idxpair])[dist], s=2, alpha=0.2)
plt.scatter(sep.to('arcsec')[dist&frat], (flux[idx]+flux[idxpair])[dist&frat], s=2, alpha=0.2)
plt.scatter(sep.to('arcsec')[msource], (flux[idx]+flux[idxpair])[msource], s=2, alpha=0.2)
plt.xlabel('separation [arcsec]')
plt.ylabel('S1+S2')




#f_nn_idx,f_nn_sep2d,_ = ac.match_coordinates_sky(c,c,nthneighbor=2)
    


### write output file

#if os.path.exists(lofarcat_file_srt):
    #os.remove(lofarcat_file_srt)
#lofarcat.write(lofarcat_file_srt)