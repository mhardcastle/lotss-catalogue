Differences for the deep fields:

-- We shouldn't need any of the remove_lgz_sources.py stuff. However
   it does change the format of the catalogue from the output of lgz_process
   so downstream code gets altered

-- We would then do zooms and blends from the lgz_process files

-- We need to process the workflow-selected blends as well as the ones
   identified by LGZ.

-- It would be best to do the whole process in one go -- from the
   workflow output to the final component and source catalogues.
   
