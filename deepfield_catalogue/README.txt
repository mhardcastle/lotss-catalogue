To create the final catalogue:

* Work in directories named for the field (en1, bootes, lockman)
* Run lgz_process/aggregate_lofgalzoo_deep.py to generate the LGZ inputs
* Run deepfield_catalogue/make_catalogue.py to ingest prefilter LGZ, blend, TZI
* Run deepfield_catalogue/merge_optical.py to merge with the optical and science-ready catalogues

To bump version number change make_catalogue, other code should then pick this up.
