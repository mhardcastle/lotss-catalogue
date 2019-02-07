New Radio Galaxy Zoo
====================

Set path somewhere sensible: set IMAGEDIR to the image dump directory

export LGZPATH=`pwd` # in the lotss-catalogue directory
export PYTHONPATH=$LGZPATH/utils:$PYTHONPATH
export IMAGEDIR=/data1/MRP1/LRGZ/images
export LOTSS_COMPONENT_CATALOGUE=components.fits

Montage needs to be on your PATH: e.g. at Herts

export PATH=/soft/Montage_v3.3/bin:$PATH

Make a FITS file for components to be used in the test: should have
roughly the same format as the input component catalogue. In
particular source names should match.

Download the image files:

python $LGZPATH/utils/download_image_files.py file.fits

Make a single image:

python $LGZPATH/lgz_create/make_overlays_legacy.py file.fits N

where N = 0 .. number of sources.

Make a set of images:

python $LGZPATH/lgz_create/make_overlays_legacy.py file.fits 0 10

Make the manifest file for upload

cat > manifest.csv <<EOF
subject_id,image_name_1,image_name_2,source_name,ra,dec,#size
EOF
cat *-manifest.txt | sort -n -k 1 -t , >> manifest.csv

