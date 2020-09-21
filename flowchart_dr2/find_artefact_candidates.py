#!/usr/bin/python
import sys
import numpy as np
from astropy.io import fits
from astropy.table import Table
import astropy.coordinates as ac

####################################




if len(sys.argv) == 1:
    print("Usage is : python find_artefact_candidates.py field_code ")
    print('E.g.: python find_artefact_candidates.py 0 ')
    sys.exit(1)

h = str(sys.argv[1])
if 'h' not in h:
    h+='h'
if h not in  ['0h','13h','n0h','n13h','s0h','s13h']:
    print('unknown field code (should be 0h or 13h)',h)
    sys.exit(1)

path = '/data2/wwilliams/projects/lofar_surveys/LoTSS-DR2-Feb2020/'
lofarcat_file_srt = path+'LoTSS_DR2_v100.srl_{h}.lr-full.presorted.fits'.format(h=h)



lofarcat = Table.read(lofarcat_file_srt)

lofarcat = Table.read(lofarcat_file_srt,format='fits')

# select sources bright and small in size
mask_bright_small = (lofarcat['Total_flux']>=5.) & (lofarcat['Maj']<=15.)

# mask the source catalogue to search for artefacts around them
source_cat = lofarcat[mask_bright_small]
source_index = source_cat['Index']
source_id = source_cat['Source_Name']

csource = ac.SkyCoord(source_cat['RA'], source_cat['DEC'], unit="deg")

# eliminate selected target sources
sea_id = set(lofarcat['Source_Name']) - set(source_id)
search_id = list(sea_id)
mask = np.in1d(lofarcat['Source_Name'],search_id)
search_cat = lofarcat[mask]
search_id = search_cat['Index']

csearch = ac.SkyCoord(search_cat['RA'], search_cat['DEC'], unit="deg")


artefact_ids = [];artefact_ra=[];artefact_dec=[]
source_ids = [];source_ra = [];source_dec = []

art_nn_idx,art_nn_sep2d,art_nn_dist3d = ac.match_coordinates_sky(csource,csearch,nthneighbor=1)


for i in art_nn_idx:
    artefact_ids.append(search_cat[i]['Source_Name'])


sep = np.array(art_nn_sep2d*3600.)
arte_index = [];seps=[];artefact_names=[]


# search for artefact candidates which have size * 1.5 of the target source and the separation 
# between target source and artefact candidate to be equal and less than 10 arcseconds.

for i,k,j in zip(artefact_ids,source_id,sep):
    if search_cat['Maj'][search_cat['Source_Name']==i]/source_cat['Maj'][source_cat['Source_Name']==k]>=1.5 and j<=10.:
        artefact_names.append(i)
        source_ids.append(k[0])    
        seps.append(j)

t=Table([artefact_names,source_ids,seps])
t.write('artefact_candidates_{h}.fits'.format(h=h),format='fits',overwrite=True)    




