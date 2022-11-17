zoom.py : deal with too-zoomed-in sources by creating new compound sources.

For bash:

export LGZPATH=/home/mjh/git/lotss-catalogue
export IMAGEDIR=/data/lofar/DR2
export PATH=/soft/Montage_v3.3/bin:$PATH

export PYTHONPATH=${LGZPATH}/utils:${LGZPATH}/dr2_catalogue:$PYTHONPATH

For (t)csh:

setenv LGZPATH /home/mjh/git/lotss-catalogue
setenv IMAGEDIR /data/lofar/DR2
setenv PATH /soft/Montage_v3.3/bin:$PATH

setenv PYTHONPATH ${LGZPATH}/utils:${LGZPATH}/dr2_catalogue:$PYTHONPATH

The Too-Zoomed-In (TZI) workflow displays an interactive
(matplotlib/aplpy) view of a source that has been marked as needing
human visual inspection by an LGZ user or by the detection of some
anomalous feature such as two or more sources sharing a single optical
ID.

The user of the TZI workflow is shown a particular source (made up of
individual components) and is asked to add or remove components and to
mark or occasionally remove an optical ID. Once this is done a 'zoom
file' is written which is a text file recording the actions to be
taken for that source. When the final catalogue is generated all of
these actions will be applied.

At the start of a run the current source is marked with a red cross
and will have one or more components indicated with green
ellipses. Components that are not part of the current source are
marked in cyan. The optical ID, if present, is marked with a magenta X.

Addition or removal of components uses clicks on the image. Clicking
on the image uses the middle button (or scroll wheel) of a three-button
mouse. This is because matplotlib grabs the other two buttons for
panning/resizing actions. What the clicks do depends on the mode
you're in, which is printed in the text window. You need to run in a
display capable of showing the text and graphics windows simultaneously.

General principles are:

-- make sure that you have a view that's capable of seeing whether the
   selected source is actually a component of a much larger
   object. Giant radio galaxies can extend for much larger than the
   field of view of the original image. Use the Z option if you need
   to expand the original cutout used.

-- select a plausible optical ID if you can, making use of the WISE
   and optical data. The latter, which is shown to RGZ users, has much
   better resolution but is less sensitive. Radio galaxy hosts are
   very often strong WISE sources but only weak optical sources
   (because they are at high redshift and dominated by an old, red
   stellar population) so the most plausible sources will have WISE
   counterparts -- unless a bright object in the field swamps the WISE
   data. Be reasonably generous if an ID is physically plausible since
   it's better for WEAVE to observe some rubbish than to miss real
   host galaxies. Don't be afraid to change an existing RGZ ID if you
   think it's necessary since these don't have the benefit of the WISE
   data.

-- To remove an optical ID just select the 'o' option but then don't
   click on anything.

-- Toggle the display of catalogued galaxies on with T if you aren't
   sure of the exact position of an optical counterpart, but don't
   select an object that you can't actually see in an image on the whole.

-- Use the FITS (f) option to fire up a ds9 window -- it's often much
   more obvious what the source is and also potentially easier to see
   other components that belong with it (since it's easier to see
   matching surface brightness levels). It may be worth setting your
   ds9 preferences to have a log scale at startup. (Run /soft/bin/ds9,
   go to Edit -> Preferences -> Menus and Buttons -> Scale -> Log,
   then click Save.)

-- Drop (d) sources that are obvious artefacts or not real sources for
   some other reason. If not obvious artefacts but they don't have any
   optical ID or  association, just leave them by saving  (s) the zoom
   file unchanged -- that keeps them in the catalogue for someone else
   to investigate  later. There  will be some  such sources  e.g. very
   high-z galaxies.
   
-- You will see a lot of sources where two ellipses, a small and a big
   one, have the same optical ID (this is one of the selection
   criteria for checking in this way). If you think the larger ellipse
   is actually representing real emission you should join them up. If
   you think the larger ellipse is an artefact then skip it and you
   will probably get the chance to drop it later.

-- When finished for a session (if you're going to be away for a
   little while) probably best to use the (q)uit option because the
   connection to the database server times out if left alone. For me
   the matplotlib window doesn't survive an x2go session close and
   restart either.

DR2 notes:

* New codes:

-- zoomprep_dr2.py updates the SQL database with the records for
   sources that need zooms and downloads any images. Should be run
   after a make_catalogue.py ... save run to create the latest input
   files.

-- zoom_dr2_sql.py uses this database to run the zooms. Relative to
   earlier versions this mostly differs by better integration with ds9
   -- a ds9 is automatically started -- and a few extra options.
   
