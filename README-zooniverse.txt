Making image sets and uploading to Zooniverse
=============================================

Set path somewhere sensible: set IMAGEDIR to the image dump directory

export LGZPATH=/home/mjh/git/lgz
export IMAGEDIR=/data/lofar/mjh/hetdex_v3

Montage needs to be on your PATH: e.g. at Herts

export PATH=/soft/Montage_v3.3/bin:$PATH

Make a source FITS file with 'RA', 'DEC', 'Source_ID' columns.

If source ID is broken, fix it -- will rewrite based on RA,DEC:

python $LGZPATH/fix_sourceid.py infile outfile 

Downloads all required files to IMAGEDIR and makes the image list as
file-list.txt

python $LGZPATH/download_image_files.py file.fits

To make a single image set:

python $LGZPATH/make_overlays.py file.fits xxxx

where xxxx is the number of the image from file.fits or file-list.txt.

To make all the images, use a job:

wc file-list.txt

(to find how many jobs you're going to be running) and then

qsub t 0-100 -v INFILE=file.fits,LGZPATH=$LGZPATH,IMAGEDIR=$IMAGEDIR $LGZPATH/lgz.qsub

To upload to Zooniverse...

export PANOPTES_PASSWORD=whatever
panoptes-subject-uploader ./manifest.csv --username mjh22 --project 2513 --workflow 1772
