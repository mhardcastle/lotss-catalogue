# Final filtering on MKSP table

import sys
from astropy.table import Table,MaskedColumn
import numpy as np

legacy_colnames=[ 'Legacy_ID', 'release', 'brickid', 'objid', 'maskbits', 'fracflux_g', 'fracflux_r', 'fracflux_z', 'type', 'Legacy_RA', 'Legacy_DEC', 'pstar', 'star', 'lupt_g', 'lupterr_g', 'lupt_r', 'lupterr_r', 'lupt_z', 'lupterr_z', 'lupt_w1', 'lupterr_w1', 'lupt_w2', 'lupterr_w2', 'lupt_w3', 'lupterr_w3', 'lupt_w4', 'lupterr_w4', 'lupt_s', 'lupterr_s', 'ANYMASK_OPT', 'gmmcomp', 'zphot', 'zphot_err', 'var.density', 'var.tr.noise', 'var.in.noise', 'flag_qual', 'Separation']

t=Table.read(sys.argv[1])

for k in legacy_colnames:
    t[k]=MaskedColumn(t[k],mask=[False]*len(t))

for i,r in enumerate(t):
    if r['Separation']>1.5:
        for k in legacy_colnames:
            t[k].mask[i]=True

zphot=t['zphot']
zphot=np.where(zphot>5,np.nan,zphot)

# redshift logic

zspec_t=np.where(t['zwarning_sdss']==0,t['zspec_sdss'],np.nan)
zsource=np.where(~np.isnan(zspec_t),"SDSS","")
zspec=np.where(np.isnan(zspec_t) & ~np.isnan(t['z_desi']),t['z_desi'],zspec_t)
zsource=np.where(np.isnan(zspec_t) & ~np.isnan(t['z_desi']),"DESI",zsource)
zspec_t=zspec
zspec=np.where(np.isnan(zspec_t) & ~np.isnan(t['z_hetdex']),t['z_hetdex'],zspec_t)
zsource=np.where(np.isnan(zspec_t) & ~np.isnan(t['z_hetdex']),"HETDEX",zsource)
zphot=np.where(t['flag_qual']==1,zphot,np.nan)
zbest_t=np.where(~np.isnan(zspec),zspec,zphot)
zsource=np.where(np.isnan(zspec) & ~np.isnan(zbest_t),"Phot",zsource)
zbest=np.where(zbest_t<0,np.nan,zbest_t)
zsource=np.where(zbest_t<0,"",zsource)
t['z_best']=zbest
t['z_source']=zsource

            
t.write(sys.argv[1].replace('.fits','-final.fits'),overwrite=True)

