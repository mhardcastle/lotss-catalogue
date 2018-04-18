import os 
from astropy.table import Table, Column, join, vstack

'''
change ID_flags
1 299803
2 966
4 3406
61 54
62 1252
63 258
311 10334
312 1694

1 -> 1
2 -> 2
4 -> 0
311 -> 31
312 -> 32
61 -> 41
62 -> 41
63 -> 42
'''


version =  '1.2'

# lofar source catalogue, gaussian catalogue and ML catalogues for each


path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'


comp_out_file = os.path.join(path,'LOFAR_HBA_T1_DR1_merge_ID_v{v:s}.comp.fits'.format(v=version))
merge_out_file = os.path.join(path,'LOFAR_HBA_T1_DR1_merge_ID_v{v:s}.fits'.format(v=version))
merge_out_full_file = merge_out_file.replace('.fits','.full.fits')
comp_out_full_file = comp_out_file.replace('.fits','.full.fits')



for f in [ comp_out_file,    comp_out_full_file, merge_out_file, merge_out_full_file]:
    print f
    cat = Table.read(f)
    
    cat['ID_flag'][cat['ID_flag']==311] = 31
    cat['ID_flag'][cat['ID_flag']==312] = 32
    cat['ID_flag'][cat['ID_flag']==61] = 41
    cat['ID_flag'][cat['ID_flag']==62] = 41
    cat['ID_flag'][cat['ID_flag']==63] = 42
    cat['ID_flag'][cat['ID_flag']==4] = 0
    
    cat.write(f, overwrite=True)