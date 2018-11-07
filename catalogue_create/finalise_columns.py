
from astropy.table import Table
import os


table_in = 'LOFAR_HBA_T1_DR1_merge_ID_optical_f_v1.2_restframe.fits'
table_out = 'LOFAR_HBA_T1_DR1_merge_ID_optical_f_v1.2a_restframe.fits'

os.system('cp -r {:s} {:s}'.format(table_in, table_out))

tt = Table.read(table_out)

if 'z_spec_source' not in tt.colnames:
    tt.rename_column('z_source', 'z_spec_source')



keeplist = ['Source_Name', 'RA', 'E_RA', 'DEC', 'E_DEC', 'Peak_flux', 'E_Peak_flux', 'Total_flux', 'E_Total_flux', 'Maj', 'E_Maj', 'Min', 'E_Min', 'DC_Maj', 'E_DC_Maj', 'DC_Min', 'E_DC_Min', 'PA', 'E_PA', 'DC_PA', 'E_DC_PA', 'Isl_rms', 'S_Code', 'Mosaic_ID', 'Masked_Fraction',
            'ID_flag', 'ID_name', 'ID_ra', 'ID_dec', 'ML_LR', 'LGZ_Size', 'LGZ_Width', 'LGZ_PA', 'LGZ_Assoc', 'LGZ_Assoc_Qual', 'LGZ_ID_Qual', 'Deblended_from', 'AllWISE', 'objID', 'gFApFlux', 'gFApFluxErr', 'gFApMag', 'gFApMagErr', 'rFApFlux', 'rFApFluxErr', 'rFApMag', 'rFApMagErr', 'iFApFlux', 'iFApFluxErr', 'iFApMag', 'iFApMagErr', 'zFApFlux', 'zFApFluxErr', 'zFApMag', 'zFApMagErr', 'yFApFlux', 'yFApFluxErr', 'yFApMag', 'yFApMagErr', 'gFKronFlux', 'gFKronFluxErr', 'rFKronFlux', 'rFKronFluxErr', 'iFKronFlux', 'iFKronFluxErr', 'zFKronFlux', 'zFKronFluxErr', 'yFKronFlux', 'yFKronFluxErr', 'w1Flux', 'w1FluxErr', 'w1Mag', 'w1MagErr', 'w2Flux', 'w2FluxErr', 'w2Mag', 'w2MagErr', 'w3Flux', 'w3FluxErr', 'w3Mag', 'w3MagErr', 'w4Flux', 'w4FluxErr', 'w4Mag', 'w4MagErr',
            'z_best', 'z_best_source', 'z_spec', 'z_spec_source', 'z1_median', 'z1_min', 'z1_max', 'z1_area', 'z2_median', 'z2_min', 'z2_max', 'z2_area',
            'specAGN', 'mqcAGN', 'XrayClass', '2RXS_ID', 'XMMSL2_ID', 'IRClass', 'EBV', 'PanSTARRS_missing', 'u_rest', 'g_rest', 'r_rest', 'i_rest', 'z_rest', 'U_rest', 'B_rest', 'V_rest', 'I_rest', 'J_rest', 'Ks_rest', 'w1_rest', 'w2_rest', 'w3_rest' ]

#tt.keep_columns(keeplist)
tt = tt[keeplist]   # keep the new order

tt.write(table_out, overwrite=True)

