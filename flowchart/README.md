# flowchart
Various scripts for the LOFAR source classification decision tree including merging catalogues from ML and LGZ outputs.

Inputs:
* LOFAR_HBA_T1_DR1_catalog_v0.99.srl.fits -- PyDBSF catalogue (with corrected source names)
* LOFAR_HBA_T1_DR1_catalog_v0.99.gaus.fits -- PyBDSF Gaussian component catalogue
* lofar_pw_pdf.fits -- ML output for sources in the PyBDSF catalogue
* lofar_gaus_pw_pdf.fits -- ML output for Gaussian in the PyBDSF Gaussian component catalogue
* HETDEX-LGZ-cat-v0.11-filtered-zooms.fits -- LGZ output catalogue
* lgz_components.txt -- LGZ output components

## Workflow:
* `find_artefact_candidates.py` -- find candidate artefacts 
* `match_2masx.py` -- match to 2MASX
* `flag_edge_sources.py` -- find the edge sources to flag as artefacts
* iterate through:
  * `lofar_source_sorter.py` -- initial classification
  * `add_gaus_info.py`  -- to calculate gaussian separations, etc and store values
  * visual confirmation of artefacts and other visually checked endpoint with `classify.py` from classify
  * `get_visual_flags.py` -- apply flags generated from several passes of classify.py on different catagories
  * `handle_m_sources.py` and `get_msource_flags.py` -- run the msource flowchart and save flags
  * `lofar_source_sorter.py` -- new classification with visual flags applied
* `write_merge_catalogues.py` -- write final catalogue, merging outputs from LGZ, ML, large galaxy strands


### `get_visual_flags.py`
Adds the following columns to the PyBDSF catalogue relating to visual classification results of subsamples selected in the flowchart:
* artefact_flag
* clustered_flag
* Lclustered_flag
* huge_faint_flag
* nhuge_faint_flag
* nhuge_2masx_flag

### `add_gaus_info.py`
more info from the Gaussians:
* Ng - the number of Gaussians making up the source
* G_LR_max - max LR from all the Gaussians
* Ng_LR_good - number of Gaussians with good LR
* Flag_G_LR_problem - True if Gaussian match is ambiguous

### `handle_m_sources.py` and `get_msource_flags.py`
Adds 
* msource[12]_flag -- saves flowchart diagnosis code
  * 0: no match
  * 1: accept ML of the source
  * 2: accept ML of the gaussian with highest ML
  * 3: deblend and accept both gaussians
  * 4: deblend workflow
  * 5: LOFAR galaxy zoo 
* MC_flag[12] -- saves flowchart box number


### `match_2masx.py`
Note: some 2MASX parameters need to be fixed with `fix_2mass_catalogue.py`. Here we add the following columns to the PyBDSF catalouge:
* 2MASX -- has a *possible* 2MASX match
* 2MASX_name -- name
* 2MASX_ra, 2MASX_dec --  position
* 2MASX_size -- 2MASX r_ext value
* 2MASX_match_large -- matches to a 2MASX galaxy with r_ext >= 60"
* 2MASX_match -- matches to a 2MASX galaxy with r_ext < 60"


### `lofar_source_sorter.py`

requires pygraphviz and graphviz to plot the flowchart

function `make_sample` can be used to generate a random subsample of a given mask, to be used in making cutout images for inspection with code from lgz_create

#### Columns added to the catalogue:
from the ML catalogues:
* cLR (has nan)
* LR - the likelihood ratio  (nan's replaced with 0)
* LR_name_wise - wise name
* LR_name_ps -  panstarrs name
* LR_ra, LR_dec -  position


flags:
* Artefact_flag - combined artefact flag from several classify routes
* ID_flag -- flag indicating identification route (modified in write_merge_catalogues)
* LGZ_flag - 1 for LGZv1 ; 2 for LGZv2 ; 20 for LGZv2_zoom; etc
* FC_flag - numeric value encodes endpoint of flowchart

Nearest neighbour (NN) details:
* NN_LR - likelihood ratio of NN
* NN_sep - distance (in arcsec) to NN
* NN_idx - catalogue index of NN
* NN4_sep - distance to the 3rd nearest neighbour
* NN5_sep - distance to the 4th nearest neighbour
* NN_Total_flux - flux of NN
* NN_Frat - flux ratio of source to NN
* NN_Maj - size of NN

Nearest not-artefact neighbour (NNC) details:
* NNC_LR - likelihood ratio of NNC
* NNC_sep - distance (in arcsec) to NNC
* NNC_idx - catalogue index of NNC
* NNC4_sep - distance to the 3rd nearest neighbour
* NNC5_sep - distance to the 4th nearest neighbour
* NNC_Total_flux - flux of NNC
* NNC_Frat - flux ratio of source to NNC
* NNC_Maj - size of NNC



## `write_merge_catalogues.py`
Merge outputs from decision tree using ML or bright galaxy optical IDs where necessary and LGZ outputs where available. LGZ components are removed from the source table and replaced by the LGZ sources (in the components table, their source names are updated). Similarly, combines multiple matches to the same bright galaxy. Removes all artefacts from both the source and components tables.

Most of the columns are dropped from the catalogues (but versions with full column sets are kept as *full.fits)

Outputs:
* LOFAR_HBA_T1_DR1_merge_ID_v*.fits -- final source list, with limited columns
* LOFAR_HBA_T1_DR1_merge_ID_v*.comp.fits -- final component list


The ID_flag in the final source list is updated/simplified given LGZ outputs
* 0 for no identification possible
* 1 for LR  identification (or lack therof) 
* 2 for bright galaxy
* 3[12] for LGZ association and identification (or lack therof) 
    * 31 for LGZ
    * 32 for LGZ zoom
* 4[12] for blend association and identification (or lack therof) 
    * 41 deblending workflow 
    * 42 LGZ deblending workflow

Some new columns are added:
* ID_name (later modified after deblending...)
* ID_ra, ID_dec - position of ID from ML or LGZ
* ML_LR - the LR value for sources with an LR match
and the LGZ columns:
* LGZ_Size
* LGZ_Width
* LGZ_PA
* LGZ_Assoc
* LGZ_Assoc_Qual
* LGZ_ID_Qual
 
The component list is the PyBDSF catalogue annotated with the final source names
* Component_Name 
* Source_Name - maps to Source_Name in the source catalogue


#### misc other scripts
* `fix_2mass_catalogue.py` -- manually fix some very wrong sizes in the 2MASX catalogue
* `get_large_galaxies.py` -- produce catalogues of the sources associated with very large optical galaxies
* `get_intersecting_sources.py` -- flag the sources that intersect/overlap
* `make_summary.py` -- combine inspection plots into single pdf file
