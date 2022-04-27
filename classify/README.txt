# Postfilter step for DR2:

* make postfilter directory and work in it

* Select all sources with LGZ_Size larger than 1 arcmin

t=Table.read('../sources-v0.4.fits')
filt=t['LGZ_Size']>60
t[filt].write('large.fits')

* Run ~/git/lotss-catalogue/utils/download_image_files_wise.py to get WISE files (will mostly exist already)

* Make sure make_overlays_postid_dr2.py points to the right directory and version of the catalogue

* Run /home/mjh/git/lotss-catalogue/lgz_create/make_overlays_postid_dr2.py

* Edit and run ~/git/lotss-catalogue/classify/populate_postfilter.py

* Edit and run ~/git/lotss-catalogue/classify/postfilter_sql.py

# Zoom postfilter

* Run ~/git/lotss-catalogue/zoom/zoomprep_dr2.py

