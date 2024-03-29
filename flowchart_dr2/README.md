# flowchart_dr2
Various scripts for the LOFAR source classification decision tree including merging catalogues from ML and LGZ outputs. Used for DR2 catalogues.

requires pygraphviz and graphviz to plot the flowchart

Inputs:
* LoTSS_DR2_v110.srl_{h}.lr-full.fits -- PyDBSF catalogue including ML output
* LoTSS_DR2_v110.gaus_{h}.lr-full.fits -- PyBDSF Gaussian component catalogue including ML output


# flags

msource_flag1 - set for multiple Gaussian sources flags msource_flag1
*  0 - no match
*  1 - accept ML of the source
*  2 - accept ML of the gaussian with highest ML
*  3 - deblend and accept both gaussians
*  4 - deblend workflow
*  5 - LOFAR galaxy zoo 
*  6 - visual check 

FC_flagN - int flowchart endpoint  for step N

ID_flag
*  -99 not set
*  -1 artefact or no priority set
*  0 - no id
*  1 - LR id or no id
*  2 - large optical galaxy
*  3 - LGZ
*  4 - visual id / prefilter (pending - should not exist after step3 after  all prefiltering is done)
*  5 - tbd
*  6 - deblend
*  7 - too zoomed in after prefilter
*  8 - Uncatalogued host after prefilter

LGZ_flag
*  1 - direct selection - large and bright
*  2 - direct selection - flowchart others
*  3 - direct selection - msource flowchart
*  4 - prefilter lgz
*  5 - prefilter deblend
*  6 - prefilter too zoomed in
 
LR_flag
*  0 - no LR id
*  1 - has LR id
*  2 - has LR gaus id

Prefilter
*  0 - Send to LGZ
*  1 - Accept ML match
*  2 - No good match
*  3 - Too zoomed in
*  4 - Artefact
*  5 - Uncatalogued host
*  6 - Blend

# steps 
Scripts need to be run in order:
1. match_2masx.py
2. set_weave_priorities.py 0
3. add_gaus_info.py 0 
4. get_nearest_neighbours.py 0 
5. get_ml_flags.py 0 
6. lofar_source_sorter_dr2.py 0 1
7. handle_m_sources_dr2.py 0 1 all
8. lofar_source_sorter_dr2.py 0 2 all
9. extract_prefilter_out.py 0  (run on lofar and copy over output... scp  lofar:/beegfs/lofar/wwilliams/lofar_surveys/DR2/LoTSS_DR2_v100.srl_0h.prefilter_outputs.fits . )
10. get_prefilter_flags.py 0
11. lofar_source_sorter_dr2.py 0 3 all

Notes: 
* 1-5 are initialisation steps that need to be run once
* 2-5 can be swapped around
* each takes as first argument the field code: 0 or 13

1. step1 is running lofar_source_sorter_dr2 for the first time
2. step2 is running it after m sources have been done
3. step3 is running after prefilter outputs have been added

### `match_2masx.py`
Produces the first presorted catalogue: LoTSS_DR2_v100.srl_{h}.lr-full.presort.hdf5

Input:
* 2MASX_dr2.fits -  2MASX table downloaded from https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=%202MASX taking only relevant columns
* LoTSS_DR2_v100.srl_{h}.lr-full.fits

Note: some 2MASX parameters need to be fixed with `fix_2mass_catalogue.py`. 

Here we add the following columns to the PyBDSF catalouge:
* 2MASX -- has a *possible* 2MASX match
* 2MASX_name -- name
* 2MASX_ra, 2MASX_dec --  position
* 2MASX_size -- 2MASX r_ext value
* 2MASX_match_large -- matches to a 2MASX galaxy with r_ext >= 60"
* 2MASX_match -- matches to a 2MASX galaxy with r_ext < 60"

### `set_weave_priorities.py`
Operates on the presorted catalogue: LoTSS_DR2_v100.srl_{h}.lr-full.presort.hdf5

For 13hr:
* WEAVE_priority1 - RA>=120 & RA<=240 & DEC>=59 & DEC<=65.5
* WEAVE_priority1a - RA>=120 & RA<=240 & DEC>=60 & DEC<=65 (was the original too small area)
* WEAVE_priority2 - RA>=110 & RA<=135 & DEC>= 26 & DEC<=41
* WEAVE_priority3 - DEC>= 39.8 & DEC<=46.05 & ~WEAVE_priority2

For 0hr:
* WEAVE_priority1 - in the legacy moc legacy-south.moc

### `add_gaus_info.py`
Operates on the presorted catalogue: LoTSS_DR2_v100.srl_{h}.lr-full.presort.hdf5

required columns in radio catlogue with LR outputs
* lr - likelihood ratio -> LR (=lr when above the threshold, =nan otherwise)
* UNWISE_OBJID - unwise match name -> LR_name_wise
* UID_L - legacy match name -> LR_name_l
* ra, dec - position of match - LR_ra, LR_dec

Need to set the log LR thresh (lLR_thresh) which is hardcoded. For DR2 this is different for the 0hr and 13hr fields, and within the 13hr field it is different above and below LR_thresh_dec (=32.375)

script adds more info from the Gaussians and calculates the Msource special cases:
* Ng - the number of Gaussians making up the source
* G_LR_max - max LR from all the Gaussians
* Ng_LR_good - number of Gaussians with good LR
* Flag_G_LR_problem - True if Gaussian match is ambiguous
* G_LR_case1 - source with lr, 1 gaus with lr diff to source
* G_LR_case2 - source with lr, 2 gaus with lr, 1 diff to source
* G_LR_case3 - 


    

### `get_nearest_neighbours.py`
Operates on the presorted catalogue: LoTSS_DR2_v100.srl_{h}.lr-full.presort.hdf5

Calculates and saves nearest neighbour distances and values.

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

### `get_ml_flags.py`
Operates on the presorted catalogue: LoTSS_DR2_v100.srl_{h}.lr-full.presort.hdf5

Adds column '0.12' as column 'ML_flag' in the presorted catalogue:
Inputs:
* LoTSS_DR2_v100.srl_{h}.lr-full.presort.hdf5
* For 0hr field: GradientBoostingClassifier_lotss_dr2_0h_pred_thresholds.fits
* For 13hr field: GradientBoostingClassifier_M2_31504_17F_TT42_B1_rd_dc_opt_tlv_high/pred_thresholds_full_13h_fixnames.fits

ML_flag - flag indicating output of ML classification for identification
* 1 for LR
* 0 for LGZ


### `lofar_source_sorter.py`
requires pygraphviz and graphviz to plot the flowchart

Contains the `Mask` class which stores the boolean mask and associated information necessary for the flowchart (parent, children, colour, label)
* function `make_sample` can be used to generate a random subsample of a given mask, to be used in making cutout images for inspection with code from lgz_create
* function `submask` will make a child of the current instance applying a new mask (join masks with AND)
* the flowchart is built by defining an initial mask, and subsequent submasks for each flowchart endpoint all stored in a list (masterlist)

Usage is : python lofar_source_sorter_dr2.py field_code step_number (weave_priority)
* weave_priority is not needed for step 1

Inputs step N+1 requires the output from step N, i.e.
* step1: LoTSS_DR2_{version}.srl_{h}.lr-full.presort.hdf5
* step2: LoTSS_DR2_{version}.srl_{h}.lr-full.sorted_step1_flux{ff:.0f}.hdf5
* step3: LoTSS_DR2_{version}.srl_{h}.lr-full.sorted_step2_flux{ff:.0f}.hdf5


Hardcoded values:
*    size_large = 15.           # in arcsec
*    separation1 = 45.          # in arcsec
*    size_huge = 25.            # in arcsec
*    fluxcut = 8.               # in mJy
*    fluxcut2 = 4.              # in mJy

Output flowchart images in pdf:
flow_13h_v100_step{step}_s15_nn45_flux4_weavepri{weavepri}.pdf

### `handle_m_sources_dr2.py`
requires pygraphviz and graphviz to plot the flowchart

Usage is : python handle_m_sources_dr2.py field_code step_number mode

field_code is one of [0, 13]
step_number should be 1 (2 is leftover from DR1)
mode is one of  [all, isol, nonisol, isol-faint, nonisol-faint] for doing all, isolated, non-isolated, faint isolated or faint non-isolated Msources.

Hardcoded values:
*   size_large = 15.   # in arcsec
*   separation1 = 45.  # in arcsec
*   size_huge = 25.    # in arcsec
*   fluxcut2 = 4.      # in mJy

Output flowchart image as pdf:
flow_13h_v100_step{step}_s15_nn45_msources_{mode}_flux4.pdf


