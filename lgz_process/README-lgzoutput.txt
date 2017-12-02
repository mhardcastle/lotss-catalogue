There are two scripts for digesting and analysing Zooniverse output.

galzoo_export.py takes a Lofar-Galaxy-Zoo classifications export file and turns it into four csv files containing lists of clicks and subject info. It is best to manually filter
the Zooniverse output to only include the correct version of the workflow before running.

aggregate_lofgalzoo.py takes these four csv files as input and analyses the click data to generate output catalogues. It includes several steps that are only required because of mismatched catalogues, and so these can be removed once a stable input catalogue exists.
