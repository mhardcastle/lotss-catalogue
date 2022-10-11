import sys
import numpy as np
from astropy.table import Table, Column
import astropy.units as u    
import astropy.coordinates as ac
from shapely.geometry import Polygon, Point
import mocpy
import pymoc


def inHETDEX(ra,dec):
    inside = np.zeros(len(ra), dtype=bool)
    
    hvert = Table.read('/home/wwilliams/data2/projects/lofar_surveys/DR1/hetdex_vertices.fits')
    hra = hvert['HETDEX_RA'][0]
    hdec = hvert['HETDEX_DEC'][0]
    
    phetdex = Polygon(np.array((hra, hdec)).T)
    for i in range(len(ra)):
        inside[i] = phetdex.contains(Point(ra[i],dec[i]))
    
    return inside



def make_polygon(ra1,ra2,dec1,dec2,delta=0.01):
    nra = int((ra2-ra1)/delta)
    ndec = int((dec2-dec1)/delta)
    ras1 = np.linspace(ra1,ra2,nra)
    decs1 = dec1*np.ones_like(ras1)
    decs2 = np.linspace(dec1,dec2,ndec)
    ras2 = ra1*np.ones_like(decs2)
    ras3 = np.linspace(ra1,ra2,nra)
    decs3 = dec2*np.ones_like(ras3)
    decs4 = np.linspace(dec1,dec2,ndec)
    ras4 = ra2*np.ones_like(decs4)
    ras,decs = np.hstack((ras1,ras2,ras3,ras4)),np.hstack((decs1,decs2,decs3,decs4))
    ras = ras[:-1]
    decs = decs[:-1]
    c = ac.SkyCoord(ras,decs,unit=u.deg,frame='icrs')
    return c

path = '/Users/w.williams/projects/lofar_surveys/DR2/'

version ='v110'

h ='13h'

    



    
m_8hr = mocpy.MOC.from_polygon_skycoord(make_polygon(110,135,26,41,delta=1),max_depth=10)

m_8hr.write('mocs/dr2_8hr.moc', format='fits',overwrite=True)   

m_13hr45 = mocpy.MOC.from_polygon_skycoord(make_polygon(130,280,40,45,delta=0.01),max_depth=10).complement()
m_13hr45.write('mocs/dr2_13h45.moc',format='fits',overwrite=True)  

m_13hr60 = mocpy.MOC.from_polygon_skycoord(make_polygon(120,240,59,65.5),max_depth=10)
m_13hr60.write('mocs/dr2_13h60.moc', format='fits',overwrite=True)   



fig = plt.figure(figsize=(15, 10))
# Define a astropy WCS easily
with World2ScreenMPL(fig, 
        fov=360 * u.deg,
        center=ac.SkyCoord(180, 0, unit='deg', frame='icrs'),
        coordsys="icrs",
        rotation=ac.Angle(0, u.degree),
        projection="AIT") as wcs:
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    # Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
    
    m_13hr60.fill(ax=ax, wcs=wcs, alpha=0.75, fill=True, color="C0", label='60 strip')
    m_13hr60.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    
    m_13hr45.fill(ax=ax, wcs=wcs, alpha=0.75, fill=True, color="C1", label='45 strip')
    m_13hr45.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    
    m_8hr.fill(ax=ax, wcs=wcs, alpha=0.75, fill=True, color="C2", label='8hr block')
    m_8hr.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    
plt.legend()
plt.xlabel('ra')
plt.ylabel('dec')
plt.title('Coverage of DR2')
plt.grid(color="black", linestyle="dotted")
plt.show()



    
h ='0h'

legacy_moc_s = mocpy.MOC.from_fits('/data2/wwilliams/projects/lofar_surveys/DR2/mocs/legacy-south.moc') 
dr2_moc = mocpy.MOC.from_fits('/data2/wwilliams/projects/lofar_surveys/DR2/mocs/dr2-moc.moc') 

 
dr2_0hr_moc = mocpy.MOC.from_fits('/data2/wwilliams/projects/lofar_surveys/DR2/mocs/dr2-0hr-moc.moc')
dr2_13hr_moc = mocpy.MOC.from_fits('/data2/wwilliams/projects/lofar_surveys/DR2/mocs/dr2-13hr-moc.moc')  

dr2_13hr_60_moc = dr2_13hr_moc.intersection(m_13hr60)
dr2_13hr_60_moc.write('mocs/dr2-13hr-weavepri1-moc.moc', format='fits',overwrite=True)   
dr2_13hr_8hr_moc = dr2_13hr_moc.intersection(m_8hr)
dr2_13hr_8hr_moc.write('mocs/dr2-13hr-weavepri2-moc.moc', format='fits',overwrite=True)   
dr2_13hr_45_moc = dr2_13hr_moc.intersection(m_13hr45)
dr2_13hr_45_moc.write('mocs/dr2-13hr-weavepri3-moc.moc', format='fits',overwrite=True)   
    
dr2_0hr_legacy_moc = dr2_0hr_moc.intersection(legacy_moc_s)
dr2_0hr_legacy_moc.write('mocs/dr2-0hr-legacy-moc.moc', format='fits',overwrite=True)   

Asky = * 4*np.pi * (180./np.pi) * (180/np.pi)
Adr2 = dr2_moc.sky_fraction * Asky

A0hr = dr2_0hr_moc.sky_fraction * Asky
A13hr = dr2_13hr_moc.sky_fraction * Asky
A13hr60 = dr2_13hr_60_moc.sky_fraction * Asky



from mocpy import MOC, World2ScreenMPL

import matplotlib.pyplot as plt
fig = plt.figure(figsize=(15, 10))
# Define a astropy WCS easily
with World2ScreenMPL(fig, 
        fov=360 * u.deg,
        center=ac.SkyCoord(0, 0, unit='deg', frame='icrs'),
        coordsys="icrs",
        rotation=ac.Angle(0, u.degree),
        projection="AIT") as wcs:
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    # Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
    
    dr2_0hr_moc.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="C0", label='DR2 0hr')
    dr2_0hr_moc.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    
    dr2_0hr_legacy_moc.fill(ax=ax, wcs=wcs, alpha=0.75, fill=True, color="C0", label='DR2 0hr - legacy')
    dr2_0hr_legacy_moc.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    
    
    dr2_13hr_moc.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="C1", label='DR2 13hr')
    dr2_13hr_moc.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    
plt.legend()
plt.xlabel('ra')
plt.ylabel('dec')
plt.title('Coverage of DR2')
plt.grid(color="black", linestyle="dotted")
plt.show()

fig = plt.figure(figsize=(15, 10))
# Define a astropy WCS easily
with World2ScreenMPL(fig, 
        fov=360 * u.deg,
        center=ac.SkyCoord(180, 0, unit='deg', frame='icrs'),
        coordsys="icrs",
        rotation=ac.Angle(0, u.degree),
        projection="AIT") as wcs:
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    # Call fill with a matplotlib axe and the `~astropy.wcs.WCS` wcs object.
    
    dr2_0hr_moc.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="C0", label='DR2 0hr')
    dr2_0hr_moc.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    
    
    dr2_13hr_moc.fill(ax=ax, wcs=wcs, alpha=0.5, fill=True, color="C1", label='DR2 13hr')
    dr2_13hr_moc.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    
    
    dr2_13hr_60_moc.fill(ax=ax, wcs=wcs, alpha=0.75, fill=True, color="C1", label='DR2 13hr 60 strip')
    dr2_13hr_60_moc.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    
    dr2_13hr_45_moc.fill(ax=ax, wcs=wcs, alpha=0.75, fill=True, color="C1", label='DR2 13hr 45 strip')
    dr2_13hr_45_moc.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    
    dr2_13hr_8hr_moc.fill(ax=ax, wcs=wcs, alpha=0.75, fill=True, color="C1", label='DR2 13hr 8hr block')
    dr2_13hr_8hr_moc.border(ax=ax, wcs=wcs, alpha=0.5, color="black")
    
plt.legend()
plt.xlabel('ra')
plt.ylabel('dec')
plt.title('Coverage of DR2')
plt.grid(color="black", linestyle="dotted")
plt.show()
