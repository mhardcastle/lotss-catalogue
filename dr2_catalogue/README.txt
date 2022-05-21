The basic workflow here is:

make_catalogue.py -- creates a version of the radio+optical ID
catalogue by merging the various input sources. By default the current
largest numbered version is overwritten, specify v... on the command
line to give a version number, specify save to write out the database
files needed for e.g. the Zoom code.

[ process_overlap.py -- multithreaded code called by make_catalogue.py to extract sources with duplicate optical IDs. ]

merge_optical.py -- merge with the base optical catalogue to assign optical IDs. Wraps around Stilts. Relatively fast.

add_lsdr8_photoz.py -- merge in the photo-zs from Ken's tables. Relatively slow.

add_physical.py -- Add physical quantities and column descriptions

check_overlap.py -- For combined release make the source table

joined_components.py -- For combined release make the components table
(note there could be some overlap here)


