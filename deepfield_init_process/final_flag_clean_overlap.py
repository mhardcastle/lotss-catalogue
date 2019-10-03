from __future__ import division
import numpy as np
import sys
import os

# For catalog matching
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import search_around_sky

from astropy.table import Table

# For MOC creation files
import pymoc.io.fits
import healpy as hp
import pymoc.util.catalog
from mocpy import MOC, WCS

# from moc_util import coords_to_hpidx, inMoc
##################################################

"""
Script to add flag_clean and FLAG_OVERLAP column to indicate bright star masking and overlapping coverage
"""


def coords_to_hpidx(ra, dec, order):
    """Convert coordinates to HEALPix indexes
    Given to list of right ascension and declination, this function computes
    the HEALPix index (in nested scheme) at each position, at the given order.
    Parameters
    ----------
    ra: array or list of floats
        The right ascensions of the sources.
    dec: array or list of floats
        The declinations of the sources.
    order: int
        HEALPix order.
    Returns
    -------
    array of int
        The HEALPix index at each position.
    """
    ra, dec = np.array(ra), np.array(dec)

    theta = 0.5 * np.pi - np.radians(dec)
    phi = np.radians(ra)
    healpix_idx = hp.ang2pix(2**order, theta, phi, nest=True)

    return healpix_idx


def inMoc(ra, dec, moc):
    """Find source position in a MOC
    Given a list of positions and a Multi Order Coverage (MOC) map, this
    function return a boolean mask with True for sources that fall inside the
    MOC and False elsewhere.
    Parameters
    ----------
    ra: array or list of floats
        The right ascensions of the sources.
    dec: array or list of floats
        The declinations of the sources.
    moc: pymoc.MOC
        The MOC read by pymoc
    Returns
    -------
    array of booleans
        The boolean mask with True for sources that fall inside the MOC.
    """
    source_healpix_cells = coords_to_hpidx(np.array(ra), np.array(dec), moc.order)

    # Array of all the HEALpix cell ids of the MOC at its maximum order.
    moc_healpix_cells = np.array(list(moc.flattened()))

    # We look for sources that are in the MOC and return the mask
    return np.in1d(source_healpix_cells, moc_healpix_cells)


######################################################################

def final_flag(field_rk,final_path,final_outpath):

    # Define the overlap bits for each filter based on the field_rk
    olap_bits = dict()

    if field_rk == "en1":
        olap_bits["i"] = 1
        olap_bits["K"] = 2
        olap_bits["sw2"] = 4
    elif field_rk == "lockman":
        olap_bits["r"] = 1
        olap_bits["sw2"] = 2
    elif field_rk == "bootes":
        olap_bits["i"] = 1
    else:
        raise ValueError("Field not defined! Choose one of 'en1', 'lockman', or 'bootes'")

    # Store paths to all of the OVERLAP MOCs
    PATH_OLAP_MOC = dict()

    PATH_OLAP_MOC["en1_i"] = "/beegfs/lofar/deepfields/ELAIS_N1_optical/final_mocs/EL_EN1_i_MOC_withadd.fits"
    PATH_OLAP_MOC["en1_K"] = "/beegfs/lofar/deepfields/ELAIS_N1_optical/final_mocs/EL_EN1_K_MOC_withadd.fits"
    PATH_OLAP_MOC["en1_sw2"] = "/beegfs/lofar/deepfields/ELAIS_N1_optical/final_mocs/EL_EN1_sw2_MOC_withadd.fits"
    PATH_OLAP_MOC["lockman_r"] = "/beegfs/lofar/deepfields/Lockman_edited_cats/optical/final_mocs/LH_r_moc_order_18_MOC_with_add.fits"
    PATH_OLAP_MOC["lockman_sw2"] = "/beegfs/lofar/deepfields/Lockman_edited_cats/optical/final_mocs/LH_sw2_moc_order_18_MOC.fits"
    PATH_OLAP_MOC["bootes_i"] = "/beegfs/lofar/deepfields/Bootes_merged_optical/final_mocs/Bootes_i_MOC.fits"

    # Store the paths to the Spitzer and optical masks in each field
    PATH_SMASK = dict()
    PATH_SMASK["en1_o"] = "/beegfs/lofar/deepfields/ELAIS_N1_optical/final_mocs/star_mask/EL_EN1_smask_asec_moc_order_18_MOC.fits"
    PATH_SMASK["en1_s"] = "/beegfs/lofar/deepfields/ELAIS_N1_optical/final_mocs/star_mask/EL_EN1_smask_chi2s_moc_order_18_MOC.fits"
    PATH_SMASK["lockman_o"] = "/beegfs/lofar/deepfields/Lockman_edited_cats/optical/final_mocs/star_mask/LH_smasklh_asec_moc_order_18_MOC.fits"
    PATH_SMASK["lockman_s"] = "/beegfs/lofar/deepfields/Lockman_edited_cats/optical/final_mocs/star_mask/LH_smasklh_chi2s_moc_order_18_MOC.fits"
    PATH_SMASK["bootes_o"] = "/beegfs/lofar/deepfields/Bootes_merged_optical/final_mocs/star_mask/Bootes_smaskb_asec_MOC.fits"
    PATH_SMASK["bootes_s"] = "/beegfs/lofar/deepfields/Bootes_merged_optical/final_mocs/star_mask/Bootes_smaskb_chi2s_MOC.fits"

    final = Table.read(final_path)

    phot_filter = list(olap_bits.keys())
    max_flag = np.sum(list(olap_bits.values()))

    # The columns which will be edited/created
    folap_col = "FLAG_OVERLAP_RADIO"
    fclean_rcol = "flag_clean_radio"

    print("Maximum flag_overlap: ", max_flag)
    print("Filters: ", phot_filter)

    # Make a temp column to store the overlap values
    final["TMP_FOLAP"] = 0

    print("Applying flag_overlap")
    for phot_band in phot_filter:
        print("For {0}, applying +{1}".format(phot_band, olap_bits[phot_band]))
        moc_o = pymoc.MOC()
        pymoc.io.fits.read_moc_fits(moc_o, PATH_OLAP_MOC["{0}_{1}".format(field_rk, phot_band)])

        in_moco = inMoc(final["RA"], final["DEC"], moc_o)
        final["TMP_FOLAP"][in_moco] += olap_bits[phot_band]
        del moc_o

    print("Total no., fraction of sources in overlap using TMP_FOLAP values: {0}, {1}".format(np.sum(final["TMP_FOLAP"] == max_flag), np.sum(final["TMP_FOLAP"] == max_flag)
    /len(final)))

    # Only update for sources that actually need to be overwritten
    #final[folap_col][final[folap_col].mask] = final["TMP_FOLAP"][final[folap_col].mask]
    # Delete the temp col
    #del final["TMP_FOLAP"]
    # rename the temp column instead...
    final['TMP_FOLAP'].name=folap_col
    
    # Now deal with the star mask for radio

    print("Now adding flag_clean for radio sources")

    # Optical mask
    MOC_MASK_PATH = PATH_SMASK[field_rk + "_o"]
    cata_moc = pymoc.MOC()
    pymoc.io.fits.read_moc_fits(cata_moc, MOC_MASK_PATH)
    print("Optical MOC area: {0} sq. deg.".format(cata_moc.area_sq_deg))

    # Filter the catalogue using the MOC
    inmask = inMoc(final["RA"], final["DEC"], cata_moc)
    print("No. of sources in optical mask MOC: {0}".format(np.sum(inmask)))

    # Spitzer mask
    MOC_MASK_PATH_chi2s = PATH_SMASK[field_rk + "_s"]
    cata_moc_chi2s = pymoc.MOC()
    pymoc.io.fits.read_moc_fits(cata_moc_chi2s, MOC_MASK_PATH_chi2s)
    print("Spitzer MOC area: {0} sq. deg.".format(cata_moc_chi2s.area_sq_deg))

    # Filter the catalogue using the MOC
    inmask_chi2s = inMoc(final["RA"], final["DEC"], cata_moc_chi2s)
    print("No. of sources in Spitzer mask MOC: {0}".format(np.sum(inmask_chi2s)))

    final[fclean_rcol] = 1
    final[fclean_rcol][inmask] = 2  	    # For optical mask only
    final[fclean_rcol][inmask_chi2s] = 3        # For Spitzer mask

    print("Overwriting the catalogue now...")
    final.write(final_outpath, overwrite=True, format='fits')
    # del final

if __name__=='__main__':
    # if called on command line guess field from working directory,
    # and allow user to specify input and output names
    dir=os.getcwd()
    field=os.path.basename(dir)
    print 'field is',field
    final_flag(field,sys.argv[1],sys.argv[2])
    
