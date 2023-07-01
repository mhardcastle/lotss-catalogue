# merge in the sizeflux measurements

from __future__ import print_function
import glob
from astropy.table import Table,join,MaskedColumn
from tqdm import tqdm
import numpy as np
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
import astropy.units as u

table=sorted(glob.glob('combined-release-v*.fits'))[-1]

print('Reading table',table)

t=Table(Table.read(table),masked=True)

sft=Table(Table.read('sizeflux/LM-size-flux.fits'),masked=True)
del(sft['RA'])
del(sft['DEC'])

tj=Table(join(t,sft,keys='Source_Name',join_type='left'),masked=True)
for k in ['Maj_LoTSS','Size_LoTSS','Total_flux_LoTSS']:
    del(tj[k])

tj['New_size'].name='LM_size'
tj['New_flux'].name='LM_flux'
tj['LM_flux']*=1000
tj['Bad_flux'].name='Bad_LM_flux'
tj['Bad_image'].name='Bad_LM_image'
tj['Bad_LM_flux'][tj['Bad_LM_flux'].mask]=False
tj['Bad_LM_image'][tj['Bad_LM_image'].mask]=False
tj['Bad_LM_flux'].fill_value=False
tj['Bad_LM_image'].fill_value=False

'''
print('Create lookup dict')
sd={}
for i,s in tqdm(enumerate(t['Source_Name']),total=len(t)):
    sd[s]=i

print('Populate new columns')
lmsize=np.ones(len(t))*np.nan
lmflux=np.ones(len(t))*np.nan
badflux=np.zeros(len(t),dtype=bool)
badimage=np.zeros(len(t),dtype=bool)

for r in tqdm(sft):
    if r['Source_Name'] in sd:
        row=sd[r['Source_Name']]
        lmsize[row]=r['New_size']
        lmflux[row]=r['New_flux']
        badflux[row]=r['Bad_flux']
        badimage[row]=r['Bad_image']
    else:
        print(r['Source_Name'],'is in sizeflux table but not catalogue')

t['LM_size']=lmsize
t['LM_flux']=1000*lmflux
t['Bad_LM_flux']=badflux
t['Bad_LM_image']=badimage
'''

t=tj

# now apply the algorithm

print('Compute best size')

filt=~t['Bad_LM_flux'] & ~t['Bad_LM_image']
filt&=t['Size_from']!="Manual"
ratio=t['LM_flux']/t['Total_flux']
filt&=(ratio>0.8) & (ratio<1.2)
filt&=(t['LAS']>30)
filt&=(t['LAS']<600)

t['LAS']=np.where(filt,t['LM_size'],t['LAS'])
t['LAS_from']=np.where(filt,"Flood-fill",t['Size_from'])

print('Compute physical sizes')

# physical sizes

z=t['z_best']

angs=cosmo.kpc_proper_per_arcmin(z)*(u.arcmin/60.0)
size=angs*t['LAS']/u.kpc
t['Size']=MaskedColumn(size,mask=t['z_best'].mask)

del(t['Size_from'])

print('Write new table')
t.write(table.replace('.fits','-LM.fits'),overwrite=True)

    
