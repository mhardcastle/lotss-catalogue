#!/usr/bin/python

import numpy as np
from astropy.io import fits
from astropy.table import Table
import astropy.coordinates as ac

####################################

path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'
lofarcat_file_srt = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.srl.fits'
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
t.write('artefact_candidates.fits',format='fits',overwrite=True)	




