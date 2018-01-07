'''
add_gaus_info.py
add some info on the gaussians to the lofar source table
- it's a bit slow to do the separation calculation for all M sources, so do it once and save it
'''
import numpy as np
from astropy.table import Table, join, Column
import astropy.coordinates as ac


if __name__=='__main__':

    ### Required INPUTS
    # lofar source catalogue, gaussian catalogue and ML catalogues for each



    path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'
    lofargcat_file = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.gaus.fixed.fits'
    lofarcat_file = path+'LOFAR_HBA_T1_DR1_catalog_v0.95_masked.srl.fixed.presort.fits'
    psmlcat_file = path+'lofar_pw.fixed.fits'
    psmlgcat_file = path+'lofar_gaus_pw.fixed.fits'



    # Gaus catalogue
    lofargcat = Table.read(lofargcat_file)
    # only relevant gaussians are in M or C sources
    lofargcat = lofargcat[lofargcat['S_Code'] != 'S']

    # Source catalogue
    lofarcat = Table.read(lofarcat_file)


    lofarcat.add_column(Column(np.ones(len(lofarcat),dtype=int), 'Ng'))
    lofarcat.add_column(Column(np.ones(len(lofarcat),dtype=float), 'G_max_sep'))

    m_S = lofarcat['S_Code'] =='S'
    minds = np.where(~m_S)[0]
    for i,sid in zip(minds, lofarcat['Source_Name'][~m_S]):
        ig = np.where(lofargcat['Source_Name']==sid)[0]
        Ng = len(ig)
        lofarcat['Ng'][i]= Ng
        
        
        gcoords = ac.SkyCoord(lofargcat['RA'][ig], lofargcat['DEC'][ig])
        _, sep, _ = gcoords.match_to_catalog_sky(gcoords, nthneighbor=Ng)
        lofarcat['G_max_sep'][i] = np.max(sep.to('arcsec').value)
        
        
    lofarcat.write(lofarcat_file, overwrite=True)
        