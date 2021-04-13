from astropy.table import Table
import numpy as np
from scipy.interpolate import interp1d

t=Table.read('optical.fits')
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
t.write('optical_sizes.fits',overwrite=True)
