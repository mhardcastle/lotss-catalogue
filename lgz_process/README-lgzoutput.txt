There are two scripts for digesting and analysing Zooniverse output.

galzoo_export.py takes a Lofar-Galaxy-Zoo classifications export file and turns it into four csv files containing lists of clicks and subject info. It is best to manually filter
the Zooniverse output to only include the correct version of the workflow before running.

aggregate_lofgalzoo.py takes these four csv files as input and analyses the click data to generate output catalogues. It includes several steps that are only required because of mismatched catalogues, and so these can be removed once a stable input catalogue exists.


galzoo_export_deep.py is a minor modification of galzoo_export.py.

aggregate_lofgalzoo_deep.py is the deep field version of
aggregate_lofgalzoo.py with code to check the field ID and match with
the appropriate catalogue. It removes the catalogue mismatch code of
aggregate_lofgalzoo.py.

galzoo_export_dr2.py is the DR2 version of galzoo_export.py and
includes some additional code to flag sources that may need to be
'rescued' due to slight mismatches between L/RGZ input catalogues and
the DR2 data.

aggregate_lofgalzoo_dr2.py is the DR2 version of
aggregate_lofgalzoo_deep.py designed to work on the LGZ/UNWISE catalogues.

aggregate_lofgalzoo_dr2_hp_mp.py is a version of this greatly sped up
by using Healpix decomposition of the optical catalogue and
multiprocessing where possible.

preprocess_optical.py must be run to decompose the optical catalogue
into healpixes before using aggregate_lofgalzoo_dr2_hp_mp.py .

Note that for DR2 the mergeGZcats.py script needs to be run for the
Fall field to allow the old LGZ classifications to be included.

Complete DR2 process is:
------------------------

filter_classifications_spring.py or fall.py -- split out the relevant classifications from the RGZ output.

galzoo_export_dr2.py Spring.csv or Fall.csv

aggregate_lofgalzoo_dr2_hp_mp.py

If necessary move Fall processed data and run mergeGZcats.py
