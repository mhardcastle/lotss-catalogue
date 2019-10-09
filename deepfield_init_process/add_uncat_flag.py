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

# field = "EN1"

null_num = -999

BLEND_CAT = dict()
BLEND_CAT["EN1"] = "/disk3/rohitk/catalogue_processing/uncat_hosts/all_deblend_pos/EN1_all_deblend_opt_click_positions.fits"
BLEND_CAT["Bootes"] = "/disk3/rohitk/catalogue_processing/uncat_hosts/all_deblend_pos/Bootes_all_deblend_opt_click_positions.fits"
BLEND_CAT["Lockman"] = "/disk3/rohitk/catalogue_processing/uncat_hosts/all_deblend_pos/Lockman_all_deblend_opt_click_positions.fits"


OPT_CAT = dict()
OPT_CAT["EN1"] = "/disk3/rohitk/ELAIS_opt_swarped/dual_analysis/combine/MASTER_catalogue/EN1_MASTER_opt_spitzer_merged_cedit_apcorr.fits"
OPT_CAT["Bootes"] = "/disk3/rohitk/Bootes/Bootes_optical/MASTER_catalogue/Bootes_MASTER_opt_spitzer_merged_forML.fits"
OPT_CAT["Lockman"] = "/disk3/rohitk/Lockman_Hole/LH_optical/combine/MASTER_catalogue/LH_MASTER_opt_spitzer_merged_forLGZ.fits"

# Add the raw optical and Spitzer catalogues here
OPT_CAT["EN1_mir"] = "/disk3/rohitk/ELAIS_opt_swarped/dual_sex/chi2salldet/combined/ccat_chi2sall_det_20phot_cedit_apcorr.fits"
OPT_CAT["Bootes_mir"] = "/disk3/rohitk/Bootes/Bootes_optical/combined/ccat_Bootes_ch2_all_aper.fits"
OPT_CAT["Lockman_mir"] = "/disk3/rohitk/Lockman_Hole/LH_optical/new_combined/ccat_chi2sall_det_20phot_cedit_apcorr.fits"

MERGED_CAT = dict()
MERGED_CAT["EN1"] = "/disk3/rohitk/ELAIS_opt_swarped/dual_analysis/combine/MASTER_catalogue/EN1_MASTER_opt_spitzer_merged_cedit_apcorr.fits"
MERGED_CAT["Bootes"] = "/disk3/rohitk/Bootes/Bootes_optical/MASTER_catalogue/Bootes_MASTER_opt_spitzer_merged.fits"
MERGED_CAT["Lockman"] = "/disk3/rohitk/Lockman_Hole/LH_optical/combine/MASTER_catalogue/LH_MASTER_opt_spitzer_merged_cedit_apcorr.fits"

master_ra = "ALPHA_J2000"
master_dec = "DELTA_J2000"

for field in list(BLEND_CAT.keys()):

    print("##### {0} #####".format(field))
    blend_cat = Table.read(BLEND_CAT[field])
    print("No. of sources in blend cat: ", len(blend_cat))
    master = Table.read(OPT_CAT[field])

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

    # Select the combination of good flags first
    if field == "EN1":
        good_opt_sources = (master["FLAG_OVERLAP"] == 7)  # & (master["flag_clean"] != 3)
    elif field == "Bootes":
        good_opt_sources = (master["FLAG_OVERLAP"] == 1) & (master["FLAG_DEEP"] != 0)  # & (master["flag_clean"] != 3)
    elif field == "Lockman":
        good_opt_sources = (master["FLAG_OVERLAP"] == 3)  # & (master["flag_clean"] != 3)

    blend_coords = SkyCoord(blend_cat["RA"], blend_cat["DEC"], unit='deg', frame='icrs')
    master_coords = SkyCoord(master[master_ra][good_opt_sources], master[master_dec][good_opt_sources], unit='deg', frame='icrs')

    # Find the NN match to each of the click positions
    ind_m, sep2d, _ = match_coordinates_sky(blend_coords, master_coords, nthneighbor=1)

    # Get the bool array into blend_cat of uncatalogues hosts
    uncat_host_bool = sep2d.arcsec > sep_cut

    print("No. of uncatalogued hosts in {0}: {1}".format(field, np.sum(uncat_host_bool)))

    # Overwrite the catalogue with this flag such that FLAG_UNCAT = 1 : sources that are uncatalogued
    blend_cat["FLAG_UNCAT"] = 0
    blend_cat["FLAG_UNCAT"][uncat_host_bool] = 1

    print("Overwriting the blend catalogue")
    blend_cat.write(BLEND_CAT[field], overwrite=True)
