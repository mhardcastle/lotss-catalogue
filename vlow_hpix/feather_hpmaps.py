import healpy as hp
import numpy as np
from copy import copy

NSIDE=8192

lofar=hp.read_map(f'/data/lofar/mjh/temperature_hp-{NSIDE}.fits')
lofar_zero=np.where(np.isnan(lofar),0,lofar)
#lowres_lofar=hp.pixelfunc.ud_grade(lofar,512)

haslam=hp.read_map('/data/lofar/mjh/haslam408_dsds_Remazeilles2014.fits')
haslam*=(408/144)**2.7
t_limit=8000
haslam[haslam>t_limit]=t_limit # reduce the effect of point sources not present in the LOFAR data

#blanked_haslam=np.where(lowres_lofar>0,haslam,0)
haslam_alm=hp.map2alm(haslam,use_pixel_weights=True)

l_copy_max=300 # 180
lmax=NSIDE*3-1
lmax_haslam=512*3-1

lofar_alm=hp.map2alm(lofar_zero,use_pixel_weights=True)
new_alm=copy(lofar_alm)

for m in range(l_copy_max):
    indices_old=hp.Alm.getidx(lmax_haslam, np.arange(m, l_copy_max), m)
    indices_new=hp.Alm.getidx(lmax,np.arange(m, l_copy_max), m)
    new_alm[indices_new]=haslam_alm[indices_old]

new_map=hp.alm2map(new_alm,nside=NSIDE)
new_map[np.isnan(lofar)]=np.nan

hp.fitsfunc.write_map(f'/data/lofar/mjh/merged-weights-level-pixweight-{NSIDE}.fits',new_map,overwrite=True)

# Now do e.g.

#java -Xmx800g -jar AladinBeta.jar -hipsgen in=/data/lofar/mjh/merged-weights-level-pixweight.fits out=hips-out-NSIDE id=LoTSS_vlow_feather

