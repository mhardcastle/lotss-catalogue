import numpy as np
from matplotlib import pyplot as plt
import glob

# For catalog matching
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import search_around_sky

from astropy.table import Table, vstack, MaskedColumn
import os
##################################################

"""
Script to read the click positions from file, cut at some seaparation and then use source
"""

field = "Bootes"

null_num = -999

BLEND_CAT = dict()
BLEND_CAT["EN1"] = "../uncat_hosts/all_deblend_pos/EN1_all_deblend_opt_click_positions.fits"
BLEND_CAT["Bootes"] = "../uncat_hosts/all_deblend_pos/Bootes_all_deblend_opt_click_positions.fits"
BLEND_CAT["Lockman"] = "../uncat_hosts/all_deblend_pos/Lockman_all_deblend_opt_click_positions.fits"


# OPT_CAT = dict()
# OPT_CAT["EN1"] = "/disk3/rohitk/ELAIS_opt_swarped/dual_analysis/combine/MASTER_catalogue/EN1_MASTER_opt_spitzer_merged_cedit_apcorr.fits"
# OPT_CAT["Bootes"] = "/disk3/rohitk/Bootes/Bootes_optical/MASTER_catalogue/Bootes_MASTER_opt_spitzer_merged_forML.fits"
# OPT_CAT["Lockman"] = "/disk3/rohitk/Lockman_Hole/LH_optical/combine/MASTER_catalogue/LH_MASTER_opt_spitzer_merged_forLGZ.fits"

# Add the raw optical and Spitzer catalogues here
OPT_CAT = dict()
OPT_CAT["EN1_mir"] = "/disk3/rohitk/ELAIS_opt_swarped/dual_sex/chi2salldet/combined/ccat_chi2sall_det_20phot_cedit_apcorr.fits"
OPT_CAT["Bootes_mir"] = "/disk3/rohitk/Bootes/Bootes_optical/combined/ccat_Bootes_ch2_all_aper_cedit_apcorr.fits"
OPT_CAT["Lockman_mir"] = "/disk3/rohitk/Lockman_Hole/LH_optical/new_combined/ccat_chi2sall_det_20phot_cedit_apcorr.fits"

OPT_CAT["EN1"] = "/disk3/rohitk/ELAIS_opt_swarped/dual_analysis/combine/MASTER_catalogue/EN1_MASTER_opt_spitzer_merged_forLGZ.fits"
OPT_CAT["Bootes"] = "/disk3/rohitk/Bootes/Bootes_optical/MASTER_catalogue/Bootes_MASTER_opt_spitzer_merged_forML.fits"
OPT_CAT["Lockman"] = "/disk3/rohitk/Lockman_Hole/LH_optical/combine/MASTER_catalogue/LH_MASTER_opt_spitzer_merged_forLGZ.fits"

#################################################################################

master_ra = "ALPHA_J2000"
master_dec = "DELTA_J2000"

# These cuts are defined based on the plot generated from find_uncat_sep.py
if field == "EN1":
    sep_cut = 1.0
    scut_raw = 1.0
elif field == "Bootes":
    sep_cut = 1.0
    scut_raw = 1.5
elif field == "Lockman":
    sep_cut = 1.0
    scut_raw = 1.0

if field == "EN1" or field == "Lockman":
    numcol = "NUMBER"
    numcols = "NUMBER_SPITZER"
elif field == "Bootes":
    numcol = "ID"
    numcols = "ID_SPITZER"
    numcolo = "ID_OPTICAL"

blend_cat = Table.read(BLEND_CAT[field])

blend_coords = SkyCoord(blend_cat["RA"], blend_cat["DEC"], unit='deg', frame='icrs')
# Get the bool array into blend_cat of uncatalogues hosts
uncat_host_bool = blend_cat["FLAG_UNCAT"] == 1

print("No. of uncatalogued hosts in {0}: {1}".format(field, np.sum(uncat_host_bool)))

# Delete variables to save space
# del master, master_coords, good_opt_sources, ind_m, sep2d

#############################################################################################
# ***** Task 2: Cross-match the uncat_host_bool sources to the raw Spitzer detected
print("***** Task 2 *****")

# Load in the RAW Spitzer detected catalogue
cata_mir = Table.read(OPT_CAT[field + "_mir"])

# Select the combination of good flags first
if field == "EN1":
    good_opt_sources = (cata_mir["FLAG_OVERLAP"] == 7)  # & (cata_mir["flag_clean"] != 3)
elif field == "Bootes":
    good_opt_sources = (cata_mir["FLAG_OVERLAP"] == 1) & (cata_mir["FLAG_DEEP"] != 0)  # & (cata_mir["flag_clean"] != 3)
elif field == "Lockman":
    good_opt_sources = (cata_mir["FLAG_OVERLAP"] == 3)  # & (cata_mir["flag_clean"] != 3)

mir_coords = SkyCoord(cata_mir[master_ra][good_opt_sources], cata_mir[master_dec][good_opt_sources], unit='deg', frame='icrs')

# Find the NN match to each of the click positions
ind_m, sep2d, _ = match_coordinates_sky(blend_coords[uncat_host_bool], mir_coords, nthneighbor=1)

# Select the separation to use to say that we have the Spitzer detected source
ind_mir = ind_m[sep2d.arcsec <= scut_raw]
print("No. of uncatalogued hosts in the *raw* Spitzer catalogue: {0}".format(len(ind_mir)))

# For these indices into the raw catalogue, check if the NUMBER_SPITZER already exists in the final catalogue
master = Table.read(OPT_CAT[field])

cata_mir_matches = cata_mir[ind_mir]
in_fincat = np.isin(cata_mir[numcols][ind_mir], master[numcols])
print("No. of uncatalogued hosts that are ALREADY in the FINAL catalogue: {0}".format(np.sum(in_fincat)))

#############################################################################################
# ***** Task 3: Now add the missed sources from raw Spitzer catalogue (not in in_fincat) and write positions of sources 
print("***** Task 3 *****")

# First convert the NUMBER_OPTICAL to a masked column
masked_num_opt = MaskedColumn(np.arange(len(cata_mir_matches)), name=numcolo, mask=np.ones(len(cata_mir_matches), dtype=bool))
cata_mir_matches[numcolo] = masked_num_opt

# Add a "NUMBER" column with -999s for these sources
cata_mir_matches[numcol] = null_num

# Get the max NUMBER value
orig_max_number = np.max(master[numcol])

"""
# So, add and stack raw Spitzer catalogue sources that are not in_fincat to the master catalogue
new_master = vstack([master, cata_mir_matches[~in_fincat]])
# Re-write the NUMBER column such that only the NUMBER values for the new sources are changed
new_master["NUMBER"][orig_max_number+1:] = np.arange(orig_max_number+1, orig_max_number + np.sum(new_master["NUMBER"] == null_num)+1)
"""

print("Writing the new subset of catalogue which were detected in the raw Spitzer catalogue!")
# new_master.write(OPT_CAT[field][:-5] + "_add_uncat.fits")
OUT_APPEND = "in_raw_Spitzer_cat"
if not os.path.exists(OUT_APPEND):
    os.makedirs(OUT_APPEND)

cata_mir_matches[~in_fincat].write("{0}/{1}_sources_raw_spitzer_to_append.fits".format(OUT_APPEND, field), overwrite=True)

# Also write the positions of sources which need forced photometry
uncat_cata = blend_cat[uncat_host_bool]
uncat_nomatch = uncat_cata[sep2d.arcsec > scut_raw]

print("No. of sources needing forced photometry: {0}".format(len(uncat_nomatch)))

OPATH = "forced_phot_positions"
if not os.path.exists(OPATH):
    os.makedirs(OPATH)

# Write the positions of these to file
uncat_nomatch.write("{0}/{1}_forced_phot_pos.txt".format(OPATH, field), format='ascii', overwrite=True)
