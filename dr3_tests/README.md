At push time, these are all the last-but-one version of these codes thanks to git rebase incompetence. To be fixed later.

= Catalogue prep =

`find_missing_hpmap_hpm.py`: Run mos2hpmap_hpm
`find_missing_hpmap.py`: Run mos2hpmap
`mos2hpmap_hpm.py`: Mosaic to HP map script
`mos2hpmap.py`: Old mosaic to HP map script for flux/rms, old mosaic style

= 6C and NVSS =

`crossmatch_6C.py`: crossmatch with 6C
`crossmatch_NVSS.py`: crossmatch with NVSS
`make_bright_nn_6C.py`: Make isolated source table for 6C
`make_bright_nn_NVSS.py`: Make isolated source table for NVSS
`make_bright_nn.py`: Make isolated source table for LoTSS

= Plotting =

`plothpmap_flux.py`: Aitoff mapping, fluxes
`plothpmap.py`: Aitoff mapping
`plot-zea-wcs.py`: ZEA plotting, Tim's version
`plot-zea-wcs-ratio-hpmap.py`: ZEA plotting, use healpix map
`plot-zea-wcs-ratio.py`: ZEA plotting of ratios from healpix map.

= Source counts =

`median_stack_ncs.py`: Median stack the sourcecounts
`mergehpmap_hpm.py`: Merge HP maps for new style mosaics
`mergehpmap.py`: Make the flux, rms and area maps from hptables.
`numbercounts.py`: Make source counts.
`sc_fitmedian_mp.py`: multiprocessing source counts script. Working version.
`sc_fitmedian.py`: old version of source counts scripts.
`simulated_numbercounts.py`: Make simulated source counts, working version.

= Misc =

`make_hips_dir.py`: Make HIPS directories for new-style mosaics
