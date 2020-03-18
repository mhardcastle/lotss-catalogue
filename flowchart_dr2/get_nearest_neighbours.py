
import numpy as np
from astropy.table import Table, join, Column
import astropy.coordinates as ac
import astropy.units as u


if __name__=='__main__':

    ### Required INPUTS
    # lofar source catalogue

    path = '/data2/wwilliams/projects/lofar_surveys/LoTSS-DR2-Feb2020/'
    lofarcat_file = path+'LoTSS_DR2_rolling.srl_0h.lr.presort.fits'


    # Source catalogue
    lofarcat = Table.read(lofarcat_file)


    # ## nearest neighbour separation

    # get nearest neighbour for all sources
    c = ac.SkyCoord(lofarcat['RA'], lofarcat['DEC'], unit="deg")
    f_nn_idx,f_nn_sep2d,_ = ac.match_coordinates_sky(c,c,nthneighbor=2)
    

    #f_nn3_idx,f_nn3_sep2d,_ = ac.match_coordinates_sky(c,c,nthneighbor=3)
    f_nn4_idx,f_nn4_sep2d,_ = ac.match_coordinates_sky(c,c,nthneighbor=4)
    f_nn5_idx,f_nn5_sep2d,_ = ac.match_coordinates_sky(c,c,nthneighbor=5)
    #f_nn6_idx,f_nn6_sep2d,_ = ac.match_coordinates_sky(c,c,nthneighbor=6)

    if 'NN_sep' not in lofarcat.colnames:
        lofarcat.add_column(Column(lofarcat['LR'][f_nn_idx], 'NN_LR'))
        lofarcat.add_column(Column(f_nn_sep2d.to(u.arcsec).value, 'NN_sep'))
        lofarcat.add_column(Column(f_nn_idx, 'NN_idx'))
        lofarcat.add_column(Column(f_nn5_sep2d.to(u.arcsec).value, 'NN5_sep'))
        lofarcat.add_column(Column(f_nn4_sep2d.to(u.arcsec).value, 'NN4_sep'))
        lofarcat.add_column(Column(lofarcat['Total_flux'][f_nn_idx], 'NN_Total_flux'))
        lofarcat.add_column(Column(lofarcat['Total_flux']/lofarcat['NN_Total_flux'], 'NN_Frat'))
        lofarcat.add_column(Column(lofarcat['Maj'][f_nn_idx], 'NN_Maj'))
        
        
    lofarcat.write(lofarcat_file, overwrite=True)
        
