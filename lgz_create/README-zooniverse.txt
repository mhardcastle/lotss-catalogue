Making image sets and uploading to Zooniverse
=============================================

Set path somewhere sensible: set IMAGEDIR to the image dump directory

export LGZPATH=/home/mjh/git/lotss-catalogue
export IMAGEDIR=/data/lofar/mjh/hetdex_v4

Montage needs to be on your PATH: e.g. at Herts

export PATH=/soft/Montage_v3.3/bin:$PATH

Make a source FITS file with 'RA', 'DEC', 'Source_Name' columns.

If source ID is broken, fix it -- will rewrite based on RA,DEC:

python $LGZPATH/utils/fix_sourceid.py infile outfile 

Downloads all required files to IMAGEDIR and makes the image list as
file-list.txt

python $LGZPATH/utils/download_image_files.py file.fits

To make a single image set:

python $LGZPATH/lgz_create/make_overlays.py file.fits xxxx

where xxxx is the number of the image from file.fits or file-list.txt .

To make all the images when you have lots of them, use a job:

wc file-list.txt

(to find how many jobs you're going to be running) and then

qsub -t 0-100 -v INFILE=file.fits,LGZPATH=$LGZPATH,IMAGEDIR=$IMAGEDIR $LGZPATH/lgz_create/lgz_v2.qsub

To upload to Zooniverse...

cat > manifest.csv <<EOF
subject_id,image_name_1,image_name_2,image_name_3,image_name_4,source_name,ra,dec,#size
EOF
cat *-manifest.txt | sort -n -k 1 -t , >> manifest.csv

export PANOPTES_PASSWORD=whatever
NODE_ENV=production panoptes-subject-uploader ./manifest.csv --username mjh22 --project 2513 --workflow 5581

To visualize catalogued sources:

qsub -t 1-146 -v INFILE=rql.fits,LGZPATH=$LGZPATH,IMAGEDIR=$IMAGEDIR $LGZPATH/lgz_create/visualize.qsub

or, better still, use the visualize_sample.sh script which will do all
the work for you.

visualize_catalogued.py filename.fits xxxx yyyy works like
make_overlays, but for visualizing sources from the
catalogue. Optionally yyyy is the number of the last image. If yyyy is
set, the code loops over all numbers between xxxx and yyyy.
Some extra environment variables control this, including PLOTPOS and NO_OVERLAY.

cat > manifest.csv <<EOF
subject_id,image_name_1,image_name_2,source_name,ra,dec,#size
EOF
cat *-manifest.txt | sort -n -k 1 -t , >> manifest.csv
