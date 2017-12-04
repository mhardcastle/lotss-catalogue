# lofar-sources
LOFAR source classification

requires pygraphviz and graphviz to plot the flowchart


function `make_sample` can be used to generate a random subsample of a given mask, to be used in making cutout images for inspection with code from
https://github.com/mhardcastle/lgz

## Columns added to the catalogue:
* Ng - the number of Gaussians making up the source
* cLR (has nan)
* LR - the likelihood ratio  (nan's replaced with 0)

### flags
* artefact - flag 1 for artefacts

### Nearest neighbour (NN) details:
* NN_LR - likelihood ratio of NN
* NN_sep - distance (in arcsec) to NN
* NN_idx - catalogue index of NN
* NN5_sep - distance to the 4th nearest neighbour
* NN_Total_flux - flux of NN
* NN_Frat - flux ratio of source to NN
* NN_Maj - size of NN


## Workflow:
* `find_artefact_candidates.py` -- find candidate artefacts 
* visual confirmation of artefacts with classify.py from lgz
* `lofar_source_sorter.py` -- 
* `get_visual_flags.py` -- apply flags generated from several passes of classify.py on different catagories
* `write_catalogues.py` -- write final catalogue, merging outputs from LGZ, ML, large galaxy strands

#### misc other scripts
* `fix_2mass_catalogue.py` -- manually fix some very wrong sizes in the 2MASX catalogue
* `get_large_galaxies.py` -- produce catalogues of the sources associated with very large optical galaxies
* `get_intersecting_sources.py` -- flag the sources that intersect/overlap
* `make_summary.py` -- combine inspection plots into single pdf file
