import numpy as np
from astropy.table import Table

"""
Script to make a "lite" version of the optical catalogue before merging with the final radio catalogue
"""

field_rk = "en1"

# The NUMBER or ID column to choose depends on the field_rk chosen
if field_rk == "en1" or field_rk == "lockman":
    cols_to_keep = ["NUMBER", "NUMBER_OPTICAL", "NUMBER_SPITZER", "X_IMAGE", "Y_IMAGE"]
elif field_rk == "bootes":
    # Bootes doesn't have X_IMAGE or Y_IMAGE - but this is fine - we aren't releasing the optical mosaics for this field anyway
    cols_to_keep = ["ID", "ID_OPTICAL", "ID_SPITZER", "FLAG_DEEP"]
else:
    raise ValueError("Field not defined! Choose one of 'en1', 'lockman', or 'bootes'")

# The full *adduncat.fits optical catalogue
OPTCAT_PATH = dict()
OPTCAT_PATH["en1"] = "/beegfs/lofar/deepfields/ELAIS_N1_optical/catalogues/correct_merging/add_uncat/EN1_MASTER_opt_spitzer_merged_cedit_apcorr_adduncat.fits"
OPTCAT_PATH["lockman"] = "/beegfs/lofar/deepfields/Lockman_edited_cats/optical/add_uncat/LH_MASTER_opt_spitzer_merged_cedit_apcorr_adduncat.fits"
OPTCAT_PATH["bootes"] = "/beegfs/lofar/deepfields/Bootes_merged_optical/add_uncat/Bootes_MASTER_opt_spitzer_merged_adduncat.fits"

print("##### {0} #####".format(field_rk))

# Columns to keep
cols_to_keep.extend(["ALPHA_J2000", "DELTA_J2000", "flag_clean", "FLAG_OVERLAP"])

print("Columns that will be kept for mergning: ")
print(cols_to_keep)

master_full = Table.read(OPTCAT_PATH[field_rk])

# Remove unnecessary columns
master_outpath = OPTCAT_PATH[field_rk][:-5] + "_lite.fits"
master_full.keep_columns(cols_to_keep)

print("Writing the lite catalogue to: {0}".format(master_outpath))
master_full.write(master_outpath, format='fits', overwrite=True)
