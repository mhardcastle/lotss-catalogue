# Preprocessing of the optical data to simplify the selection of matches

from __future__ import print_function
from astropy_healpix import HEALPix
from astropy import units as u
import numpy as np
from astropy.table import Table
from scipy.interpolate import interp1d
from multiprocessing import Pool
from tqdm import tqdm

hp = HEALPix(nside=256)

print('HP pixels are',hp.npix)
print('HP pixel area is',60*60*hp.pixel_area.value*(180/np.pi)**2,'sq arcmin')

t=Table.read('source_lr.fits')

lhp=hp.lonlat_to_healpix(t['RA'],t['DEC'])

lpix=sorted(list(set(lhp)))

print('There are',len(lpix),'LOFAR pixels')

print('Loading optical catalogue')

optcat='optical.fits'
opt=Table.read(optcat)
idname='ID'
#opt['ID']=np.where(opt['UNWISE_OBJID']!="N/A",opt['UNWISE_OBJID'],opt['UID_L'])
rast='RA'
decst='DEC'

print('Instrumenting it with sizes')

# instrument optical table with sizes
t=opt
ts=Table.read('sz_mag.dat',format='ascii.commented_header')

rc=3.0 # default circular radius
column='sz_90'
rmax=np.max(ts[column]) # largest value for bright objects
rmax=np.sqrt(rc**2.0+rmax**2.0)
magmin=np.min(ts['mag_r'])
magmax=np.max(ts['mag_r'])

ifn=interp1d(ts['mag_r'],ts[column],fill_value='extrapolate')

rf=ifn(t['MAG_R'])

rf=np.sqrt(rf**2.0+rc**2.0)
rf=np.where(np.isnan(t['MAG_R']),rc,rf)
rf=np.where(t['MAG_R']<magmin,rmax,rf)
rf=np.where(t['MAG_R']>magmax,rc,rf)

t['size']=rf
opt=t

print('Computing healpixes')

ohp=hp.lonlat_to_healpix(opt[rast]*u.deg,opt[decst]*u.deg)

print('Multiprocessing to generate sub-cats')

def filt(pix):
    return (pix,opt[ohp==pix])

def write(result):
    pix,cat=result
    cat.write('optical_hpix_256/%i.fits' % pix,overwrite=True)

p=Pool(96)
subcats=[]
for result in tqdm(p.imap(filt,lpix),total=len(lpix)):
    subcats.append(result)
    
#subcats=p.map(filt,lpix)

print('Writing out sub-catalogues')

for _ in tqdm(p.imap(write,subcats),total=len(lpix)):
    pass

p.close()
print('Done')
