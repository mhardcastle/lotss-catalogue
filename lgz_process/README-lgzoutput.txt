There are two scripts for digesting and analysing Zooniverse output.

galzoo_export.py takes a Lofar-Galaxy-Zoo classifications export file and turns it into four csv files containing lists of clicks and subject info. It is best to manually filter
the Zooniverse output to only include the correct version of the workflow before running.

aggregate_lofgalzoo.py takes these four csv files as input and analyses the click data to generate output catalogues. It includes several steps that are only required because of mismatched catalogues, and so these can be removed once a stable input catalogue exists.


galzoo_export_deep.py is a minor modification of galzoo_export.py.

aggregate_lofgalzoo_deep.py is the deep field version of
aggregate_lofgalzoo.py with code to check the field ID and match with
the appropriate catalogue. It removes the catalogue mismatch code of
aggregate_lofgalzoo.py.

aggregate_lofgalzoo_dr2.py is the DR2 version of
aggregate_lofgalzoo_deep.py designed to work on the LGZ/UNWISE catalogues.
