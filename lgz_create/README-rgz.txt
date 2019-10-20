New Radio Galaxy Zoo
====================

Set path somewhere sensible: set IMAGEDIR to the image dump directory

export LGZPATH=/data1/MRP1/LRGZ/Hatfield/lotss-catalogue #`pwd` in the lotss-catalogue directory
export PYTHONPATH=$LGZPATH/utils:$PYTHONPATH
export IMAGEDIR=/data1/MRP1/LRGZ/images
export LOTSS_COMPONENT_CATALOGUE=/data1/MRP1/LRGZ/components/LOFAR_HBA_T1_DR1_merge_ID_v1.2.comp.fits

(for Herts, IMAGEDIR=/data/lofar/DR2 )

Montage needs to be on your PATH: e.g. at Herts

export PATH=/soft/Montage_v3.3/bin:$PATH

Make a FITS file for components to be used in the test: should have
roughly the same format as the input component catalogue. In
particular source names should match.

Download the image files:
(make sure that the -list.txt file is deleted)
# Edit notes: 
# I commented out lines 26 and 27 in download_image_files_legacy.py and lines 34 and 160 in make_overlays_legacy.py. The reason for this was because I was getting an EOF error in the get_wise function when running the download script. So I decided to remove all the downloads and further processing of the WISE and FIRST data. This may need to be reverted if we ever decide to use this data. 

python $LGZPATH/utils/download_image_files_legacy.py file.fits

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

Upload:

export PANOPTES_PASSWORD=whatever
NODE_ENV=production panoptes-subject-uploader ./manifest.csv --username mjh22 --project 8190 --workflow 11973

(or for the LGZ version)

NODE_ENV=production panoptes-subject-uploader ./manifest.csv --username mjh22 --project 2513 --workflow 12374
