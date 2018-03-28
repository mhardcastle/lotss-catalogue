Utility routines used by other parts of the catalogue codebase.

This directory needs to be on the PYTHONPATH.

All of these routines assume that $IMAGEDIR points to a directory with subdirectory 'mosaics'
(containing the LOFAR mosaics) and 'downloads' (containing downloaded image files).

* download_image_files.py : Given a FITS catalogue, downloads the FIRST,
          PanSTARRS and WISE images. Creates a standard-format text file used
	  by all the visualization code.

* fix_sourceid.py : old routine to fix up broken source IDs.

* get_fits.py : make a FITS cutout for a given source name or at a given RA, DEC

* image_utils.py : utilities related to LOFAR images used by other code

* make_image_list.py : as for download_image_files.py, but works on
          the assumption that all the relevant files have already been
          downloaded. Attempts to find a dump of all the maps and
          their parameters in a data structure
          $IMAGEDIR/mapslist.pickle -- if this does not exist, it will
          be created.

* overlay.py : the radio/optical overlay code used by LGZ and other
          visualization routines, along with the zoom, blend etc
          codes. Highly configurable, complicated code which needs
          more documentation than this.

* separation.py : position separation code.

* subim.py : code to flatten a radio FITS file to 2D and extract a subimage of given size from it.
