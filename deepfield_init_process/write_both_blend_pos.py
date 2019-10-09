import numpy as np
from matplotlib import pyplot as plt
import glob

# For catalog matching
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import search_around_sky

from astropy.table import Table
import os
##################################################

field_struct = "Bootes"

# Add sources from the four possible "deblend" outcomes:
"""
# The three directories contain output .txt files from:
1. deblends from workflow
2. deblends from LGZ
3. deblends from pre-filter
4. uncat hosts from pre-filter
"""
txt_files = glob.glob("/disk3/rohitk/catalogue_processing/uncat_hosts/prefilt_blends/{0}/*.txt".format(field_struct))
txt_files.extend(glob.glob("/disk3/rohitk/catalogue_processing/uncat_hosts/blends/{0}/*.txt".format(field_struct)))
txt_files.extend(glob.glob("/disk3/rohitk/catalogue_processing/uncat_hosts/lgz/{0}/*.txt".format(field_struct)))

print("## {0} ##".format(field_struct))
print("No. of sources in this catalogue: ", len(txt_files))

# dir_name = path_txt.split("/")[-2]
dir_name = field_struct

if dir_name == "ELAIS_N1":
    field = "EN1"
elif dir_name == "Bootes":
    field = "Bootes"
elif dir_name == "Lockman":
    field = "Lockman"
else:
    raise ValueError("Field not defined!")
    quit


# List to store the Source_Name for all potential blends that were sent to deblending workflow
ilt_sent_to_blend = []
all_opt_ids = []

for fname in txt_files:
    sname_fname = fname.split("/")[-1].split(".txt")[0]
    # Read in each blend output
    with open(fname, "r") as fin:
        lines = fin.read().splitlines()

    for ii, line in enumerate(lines):
        if line.startswith("## Optical IDs"):
            opt_ids = [kk.split(" ")[1:] for kk in lines[ii+1:]]
            all_opt_ids.extend(opt_ids)
            ilt_sent_to_blend.extend(len(opt_ids) * [sname_fname])

# Get the array of optical IDs
opt_id_arr = np.array(all_opt_ids).astype(float)

T = Table()
T["RA"] = opt_id_arr[:, 0]
T["RA"].unit = u.deg
T["DEC"] = opt_id_arr[:, 1]
T["DEC"].unit = u.deg
T["Source_Name"] = np.array(ilt_sent_to_blend)

if not os.path.exists("all_deblend_pos"):
    os.makedirs("all_deblend_pos")

print("Writing click positions from blends/ and prefilt_blends/ to: /disk3/rohitk/catalogue_processing/uncat_hosts/all_deblend_pos/{0}_all_deblend_opt_click_positions.fits".format(field))
T.write("/disk3/rohitk/catalogue_processing/uncat_hosts/all_deblend_pos/{0}_all_deblend_opt_click_positions.fits".format(field), format='fits')
