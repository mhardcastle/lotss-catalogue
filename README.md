# lotss-catalogue

This repository contains the scripts needed to generate the LOTSS
optically identified catalogue.

Note that with a few exceptions these contain hard-wired paths and are intended to be run by particular people on the LOFAR-UK (Herts) system. It may be possible to improve this situation later.

Please use this repository to open issues related to the cataloguing.

Directories contain README files which give more information on what scripts to run in what order.

* catalogue_create: mostly concerned with merging LGZ and flowchart outputs and adding the optical IDs

* classify: code allowing interactive user classification of sources

* lgz_create: code for making the LGZ images

* lgz_process: Judith's LGZ output processing code

* old: some code no longer or not currently used

* utils: generally useful code which needs to be on the PYTHONPATH for some other routines

* zoom: interactively deal with too-zoomed-in sources.
