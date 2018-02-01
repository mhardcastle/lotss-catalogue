Routines in this directory create the catalogue from the
LOFAR-Galaxy-zoo and maximum-likelihood catalogues.

remove_lgz_sources.py : sorts out the matching between LGZv1 outputs
and the actual source catalogue. Makes components.txt and remove.txt files.

../zoom/zoom.py : Use this on any file made from the LGZ catalogue
(but specifically intended for those with high Zoom_prob) to generate
override files in the zoom directory.

process_hostbroken.py : sort out the host_broken sources

process_lgz.py : makes a combined LGZ and too-zoomed-in catalogue from
the output of remove_lgz_source.py and the files in the zoom directory.

[Wendy's flowchart code runs next merging the LGZ catalogue and making
use of the new_components.txt file]

process_blends.py : operate on the blend files to sort out the blends

Topcat/stilts merge Wendy's catalogue with the optical ID catalogue, then

fixup_topcat_merge.py : rename and remove columns to give the final catalogue
