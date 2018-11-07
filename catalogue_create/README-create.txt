Routines in this directory create the catalogue from the
LOFAR-Galaxy-zoo and maximum-likelihood catalogues.

remove_lgz_sources.py : sorts out the matching between LGZv1/2 outputs
and the actual source catalogue. Makes components.txt and remove.txt files.

../zoom/zoom.py : Use this on any file made from the LGZ catalogue
(but specifically intended for those with high Zoom_prob) to generate
override files in the zoom directory.

process_hostbroken.py : sort out the host_broken sources

process_lgz.py : makes a combined LGZ and too-zoomed-in catalogue from
the output of remove_lgz_sources.py and the files in the zoom directory.

[Wendy's flowchart code runs next merging the LGZ catalogue and making
use of the new_components.txt file]

In blend directory:

../blend/extract_blends.py : make blends catalogue

../utils/make_image_list.py blends.fits : make image file for catalogue

../blend/blend_components.py blends.fits : use on a file of potential blended
sources to generate override files in the blends directory.

process_blend.py : operate on the blend files and output from Wendy's code to sort out the blends

In output directory:

process_2mass_tzi.py : a few 2MASS sources come through the
too-zoomed-in path. Fix up their ra, dec manually and change their
ID-flag so that they keep their names.

Topcat/stilts merge blend output catalogue with the optical ID catalogue, then

fixup_topcat_merge.py : rename and remove columns to give the final catalogue
