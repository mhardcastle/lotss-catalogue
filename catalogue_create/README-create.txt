Routines in this directory create the catalogue from the
LOFAR-Galaxy-zoo and maximum-likelihood catalogues.

remove_lgz_sources.py : sorts out the matching between LGZv1 outputs
and the actual source catalogue.

process_lgz.py : makes a combined LGZ and too-zoomed-in catalogue from
the output of remove_lgz_source.py.

[Wendy's flowchart code runs next merging the LGZ catalogue and making
use of the remove.txt file]

Topcat/stilts merge Wendy's catalogue with the optical ID catalogue, then

fixup_topcat_merge.py : rename and remove columns to give the final catalogue



