from __future__ import print_function
import sys
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
import astropy.units as u
import numpy as np
from astropy.io import fits

infile=sys.argv[1]

t=Table.read(infile)
tg=Table.read('/beegfs/lofar/mjh/rgz/gloudemans_match.fits')

# add resolved column using the Shimwell+ 21 criterion

snr=t['Total_flux']/t['E_Total_flux']
r999=0.42+(1.08/(1+(snr/96.57)**2.49))
resolved=np.log(t['Total_flux']/t['Peak_flux'])>r999
resolved|=(t['S_Code']=='Z')
t['Resolved']=resolved
asize=np.where(t['Composite_Size']>0,t['Composite_Size'],2*t['DC_Maj'])
t['LAS']=asize

zphot=t['zphot']
zphot=np.where(zphot>5,np.nan,zphot)

# redshift logic

zspec_t=np.where(t['zwarning_sdss']==0,t['zspec_sdss'],np.nan)
zsource=np.where(~np.isnan(zspec_t),"SDSS","")
zspec=np.where(np.isnan(zspec_t) & ~np.isnan(t['z_hetdex']),t['z_hetdex'],zspec_t)
zsource=np.where(np.isnan(zspec_t) & ~np.isnan(t['z_hetdex']),"HETDEX",zsource)
zphot=np.where(t['flag_qual']==1,zphot,np.nan)
zbest_t=np.where(~np.isnan(zspec),zspec,zphot)
zsource=np.where(np.isnan(zspec) & ~np.isnan(zbest_t),"Phot",zsource)
zbest=np.where(zbest_t<0,np.nan,zbest_t)
zsource=np.where(zbest_t<0,"",zsource)
t['z_best']=zbest
t['z_source']=zsource

# merge in the Gloudemans+ high-z quasar redshifts
for r in tg:
    print('Doing',r['Source_Name'])
    mask=(t['Source_Name']==r['Source_Name'])
    if np.any(mask):
        i=np.argmax(mask)
        print('Inserting in row',i)
        t[i]['z_best']=r['zsp']
        t[i]['z_source']='HZQ'

z=t['z_best']

ld=cosmo.luminosity_distance(z)
lr=4*1e-29*np.pi*t['Total_flux']*ld.to(u.m)**2.0*(1+z)**(-0.3)
angs=cosmo.kpc_proper_per_arcmin(z)*(u.arcmin/60.0)
size=angs*asize/u.kpc
t['Size']=size
t['L_144']=lr

# Add some column descriptions

lines=open('../README.md').readlines()
for l in lines:
    l=l.replace('\\','')
    if '|' not in l:
        continue
    if 'Column' in l:
        continue
    if 'Flag value' in l:
        break
    bits=l.split('|')
    if len(bits)>1:
        c=bits[1].lstrip().rstrip()
        u=bits[2].lstrip().rstrip()
        d=bits[3].lstrip().rstrip()
        if c[0]==':': continue
        t[c].description=d
        t[c].units=u

outname=infile.replace('.fits','-physical.fits')
t.write(outname,overwrite=True)

# Now work round astropy's inability to write comments into FITS tables, sigh

hdu=fits.open(outname)

nf=hdu[1].header['TFIELDS']
for i in range(1,nf+1):
    c=hdu[1].header['TTYPE%i' %i]
    try:
        hdu[1].header['TCOMM%i' %i]=t[c].description
        hdu[1].header['TUNIT%i' %i]=t[c].units
    except AttributeError:
        print(c,'has no units or description')

hdu.writeto(outname,overwrite=True)
