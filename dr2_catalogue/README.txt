The basic workflow here is:

In the main directory:

filter_classifications_XXXX.py -- start from the flowchart output and RGZ stuff and create source_lr.fits, gaussian_lr.fits

In the field directory:

Then see ../lgz-process/README-lgzoutput.txt for (RL)GZ
processing.

If necessary, remake the weights files:

user_quality_stats.py (generates consensus_frac.json)

Then rerun the (RL)GZ processing with the weighted version.

Fall field needs a special processing step to merge in the
LGZ classifications

python ~/git/lotss-catalogue/lgz_process/mergeGZcats.py

merge_prefilter.py -- get the SQL prefilter classifications in.
mv source_lr_prefilter.fits source_lr.fits

[ if necessary, remove or rename ridgeline directory ]

make_catalogue.py -- creates a version of the radio+optical ID
catalogue by merging the various input sources. By default the current
largest numbered version is overwritten: specify v... on the command
line to give a version number, specify save to write out the database
files needed for e.g. the Zoom code.

[ if necessary, run ridgeline code, remake... ]

[ process_overlap.py -- multithreaded code called by make_catalogue.py to extract sources with duplicate optical IDs. ]

merge_optical.py -- merge with the base optical catalogue to assign optical IDs. Wraps around Stilts. Relatively fast.

add_lsdr8_photoz.py -- merge in the photo-zs from Ken's tables. Relatively slow.

add_physical.py -- Add physical quantities and column descriptions

In the parent directory:

make_joined.py -- Make the joined source table (Spring and Fall) and components table.

merge_sizeflux.py -- Merge in the LM sizes and fluxes

merge_mass_catalogue.py
