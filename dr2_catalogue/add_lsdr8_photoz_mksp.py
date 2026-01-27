# K Duncan - Jun 2021
# Cross-matching LoTSS DR2 Catalogues with LS DR8 Photo-z
#
# This version works with an arbitrary input file and is to be used
# for matching to the polarized source catalogue for DR3.

from __future__ import print_function
import os
import numpy as np
import healpy
import glob
import sys

import astropy.units as u
from astropy.table import MaskedColumn, Table, join, vstack, hstack
from astropy.coordinates import SkyCoord

"""
Input files/paths
"""

infile=lofar_cat_path=sys.argv[1]
ID_RA=sys.argv[2]
ID_DEC=sys.argv[3]

t=Table.read(infile)

spec_sep = 1.5 # arcsec - matching radius for spec-zs
xray_sep = 3 # arcsec - matching radius for X-ray cross-IDs

verbose = True

# work out what hemisphere(s) we're using

critical_dec=32.375 # 32.275
    
hemispheres=[]

print('DEC range is',np.nanmin(t[ID_DEC]),np.nanmax(t[ID_DEC]))

if np.nanmax(t[ID_DEC])>=critical_dec:
    hemispheres.append('north')
if np.nanmin(t[ID_DEC])<critical_dec:
    hemispheres.append('south')

del(t)

final_tables=[]
for hemisphere in hemispheres:
    print('======= Doing hemisphere %s =======' % hemisphere)
    ancillary_data = '/beegfs/lofar/duncan/ancillary_data'

    photoz_cat_path = '/beegfs/lofar/duncan/photoz/{}_merged'.format(hemisphere)

    lofar_cat=Table.read(lofar_cat_path)
    filter_dec=lofar_cat[ID_DEC]
    # Set paths/variables given inputs
    if hemisphere == 'north':
        lofar_cat=lofar_cat[filter_dec>=critical_dec]
        release = 8001  # North = 8001, South = 8000
    else:
        lofar_cat=lofar_cat[filter_dec<critical_dec]
        release = 8000

    sdss_gal_path = '{0}/{1}'.format(ancillary_data, 'specObj-dr16.fits')
    sdss_qso_path = '{0}/{1}'.format(ancillary_data, 'DR16Q_v4.fits')
    hetdex_path = '/beegfs/lofar/mjh/rgz/hetdex_sc1_v3.2.ecsv'
    desi_path='/beegfs/lofar/mjh/rgz/zall-pix-iron.fits'
    xray_2rxs_path = '{0}/{1}'.format(ancillary_data, '2RXS_AllWISE_catalog_paper_2017May26.fits')
    xray_xmmsl2_path = '{0}/{1}'.format(ancillary_data, 'XMMSL2_AllWISE_catalog_paper_2017JUN09.fits')
    bricks_path = '{0}/{1}'.format(ancillary_data, 'survey-bricks.fits')

    """
    Merge phot-z catalogue
    """

    ra=lofar_cat[ID_RA]
    dec=lofar_cat[ID_DEC]
    
    lofar_hpx_opt = healpy.ang2pix(2**3, ra, dec, lonlat=True)

    lofar_cat['HPX']=lofar_hpx_opt

    lofar_hpx_unique = np.unique(lofar_hpx_opt)

    ls_joined_all = []

    # Loop through Healpix chunks to join photo-z information
    for hpx in lofar_hpx_unique[:]:
        if verbose: print('Healpix number is',hpx)
        lofar_subset = (lofar_hpx_opt == hpx)
        lofar_subset_cat=lofar_cat[lofar_subset]
        opt_coord=SkyCoord(lofar_subset_cat[ID_RA],lofar_subset_cat[ID_DEC], unit='deg')
        #hpx_opt = healpy.ang2pix(2**3, lofar_cat[ID_RA][with_match*lofar_subset],
        #                               lofar_cat[ID_DEC][with_match*lofar_subset],
        #                               lonlat=True)

        opt_hpx_unq = lofar_hpx_unique
        hpx_list = np.append(hpx, healpy.get_all_neighbours(2**3, hpx))
        hpx_list = hpx_list[np.in1d(hpx_list, opt_hpx_unq)] # Crop surplus neighbours

        pzpath = '{0}/hpx_{1:03}_merged.fits'
        filenames=[pzpath.format(photoz_cat_path, o) for o in hpx_list]
        if verbose:
            print('Filenames are',filenames)
        tables=[Table.read(f) for f in filenames if os.path.isfile(f)]
        if len(tables):
            print('Number of files to load are',len(tables))
            try:
                photoz = vstack(tables)
            except TypeError:
                print(tables)
                raise

            photoz['id'].name = 'Legacy_ID'
            photoz['ra'].name = 'Legacy_RA'
            photoz['dec'].name = 'Legacy_DEC'

            # Now match not on UID but on position
            
            pz_coord=SkyCoord(photoz['Legacy_RA'],photoz['Legacy_DEC'],unit='deg')
            id_pz, d2d, _ = opt_coord.match_to_catalog_sky(pz_coord)
            print('id_pz is',id_pz)
            print(type(id_pz))
            print('length of lofar table is',len(lofar_subset_cat))
            print('length of id_pz is',len(id_pz))
            match_pz = (d2d < spec_sep*u.arcsec)

            # Now we have all the ids for the lofar subset

            for k in photoz.colnames:
                column=photoz[k][id_pz]
                lofar_subset_cat[k]=column
            lofar_subset_cat['Separation']=d2d.to(u.arcsec)

            ls_joined_all.append(lofar_subset_cat)

    merged_all = vstack(ls_joined_all)
    merged_all.sort(ID_RA)

    with_match = np.where(merged_all[ID_DEC] > 0)[0]
    opt_coord = SkyCoord(merged_all[ID_RA], merged_all[ID_DEC], unit='deg')

    """
    SDSS Matching
    """

    print('SDSS galaxy matching')

    sdss_gal = Table.read(sdss_gal_path)
    sdss_gal = sdss_gal[sdss_gal['Z'] < 2.] # Limit to trustworthy non-QSO redshifts
    sdss_gal_coord = SkyCoord(sdss_gal['PLUG_RA'], sdss_gal['PLUG_DEC'], unit='deg')

    id_sdss, d2d, _ = opt_coord[with_match].match_to_catalog_sky(sdss_gal_coord)
    match_gal = (d2d < spec_sep*u.arcsec)

    sdss_z = np.ones(len(merged_all)) * np.nan
    #sdss_zqso = np.ones(len(merged_all)) * np.nan
    sdss_zwarn = np.ones(len(merged_all)) * -1
    sdss_plate = np.ones(len(merged_all)) * -1
    sdss_mjd = np.ones(len(merged_all)) * -1
    sdss_fiberid = np.ones(len(merged_all), dtype='int') * -1

    sdss_z[with_match[match_gal]] = sdss_gal['Z'][id_sdss[match_gal]]
    sdss_zwarn[with_match[match_gal]] = sdss_gal['ZWARNING'][id_sdss[match_gal]]
    sdss_plate[with_match[match_gal]] = sdss_gal['PLATE'][id_sdss[match_gal]]
    sdss_mjd[with_match[match_gal]] = sdss_gal['MJD'][id_sdss[match_gal]]
    sdss_fiberid[with_match[match_gal]] = sdss_gal['FIBERID'][id_sdss[match_gal]]

    sdss_qso = Table.read(sdss_qso_path)
    sdss_qso=sdss_qso[sdss_qso['Z']<5.0] # Remove high-z QSO per Ken advice
    sdss_qso_coord = SkyCoord(sdss_qso['RA'], sdss_qso['DEC'], unit='deg')

    id_sdss_qso, d2d, _ = opt_coord[with_match].match_to_catalog_sky(sdss_qso_coord)
    match_qso = (d2d < spec_sep*u.arcsec)

    sdss_z[with_match[match_qso]] = sdss_qso['Z'][id_sdss_qso[match_qso]]
    sdss_zwarn[with_match[match_qso]] = sdss_qso['ZWARNING'][id_sdss_qso[match_qso]]
    sdss_plate[with_match[match_qso]] = sdss_qso['PLATE'][id_sdss_qso[match_qso]]
    sdss_mjd[with_match[match_qso]] = sdss_qso['MJD'][id_sdss_qso[match_qso]]
    sdss_fiberid[with_match[match_qso]] = sdss_qso['FIBERID'][id_sdss_qso[match_qso]]

    merged_all.add_column(MaskedColumn(name= 'zspec_sdss', data=sdss_z,
                                       mask=(sdss_z == np.nan), dtype=sdss_gal['Z'].dtype))
    merged_all.add_column(MaskedColumn(name= 'zwarning_sdss', data=sdss_zwarn,
                                       mask=(sdss_zwarn == -1), dtype=sdss_gal['ZWARNING'].dtype))
    merged_all.add_column(MaskedColumn(name= 'plate_sdss', data=sdss_plate,
                                       mask=(sdss_plate == -1), dtype=sdss_gal['PLATE'].dtype))
    merged_all.add_column(MaskedColumn(name= 'mjd_sdss', data=sdss_mjd,
                                       mask=(sdss_mjd == -1), dtype=sdss_gal['MJD'].dtype))
    merged_all.add_column(MaskedColumn(name= 'fiberid_sdss', data=sdss_fiberid,
                                       mask=(sdss_fiberid == -1), dtype=sdss_gal['FIBERID'].dtype))

    """
    Merge HETDEX
    """
    print('Merging HETDEX spectra')
    
    hetdex=Table.read(hetdex_path,format='ascii.ecsv')
    hetdex=hetdex[hetdex['z_hetdex']>0]
    hetdex_coord = SkyCoord(hetdex['RA'].data, hetdex['DEC'].data, unit='deg')

    id_hetdex, d2d, _ = opt_coord[with_match].match_to_catalog_sky(hetdex_coord)
    match_gal = (d2d < spec_sep*u.arcsec)

    hetdex_z = np.ones(len(merged_all)) * np.nan
    hetdex_zconf = np.ones(len(merged_all)) * np.nan
    hetdex_sourceid = np.ones(len(merged_all), dtype='int') * -1
    
    hetdex_z[with_match[match_gal]] = hetdex['z_hetdex'][id_hetdex[match_gal]]
    hetdex_zconf[with_match[match_gal]] = hetdex['z_hetdex_conf'][id_hetdex[match_gal]]
    hetdex_sourceid[with_match[match_gal]] = hetdex['source_id'][id_hetdex[match_gal]]

    merged_all.add_column(MaskedColumn(name='z_hetdex', data=hetdex_z,mask=(hetdex_z == np.nan), dtype=hetdex['z_hetdex'].dtype))
    merged_all.add_column(MaskedColumn(name='z_hetdex_conf', data=hetdex_zconf,mask=(hetdex_zconf == np.nan), dtype=hetdex['z_hetdex_conf'].dtype))
    merged_all.add_column(MaskedColumn(name='hetdex_sourceid', data=hetdex_sourceid,mask=(hetdex_sourceid == -1), dtype=hetdex['source_id'].dtype))
                                                                 

    """
    Merge DESI
    """
    print('Merging DESI spectra')
    
    desi=Table.read(desi_path)
    filt=desi['ZCAT_PRIMARY'] & (desi['ZWARN']==0)
    filt&=(desi['Z']<1.5) | (desi['SPECTYPE']!='GALAXY')
    desi=desi[filt]
    desi_coord = SkyCoord(desi['TARGET_RA'].data, desi['TARGET_DEC'].data, unit='deg')

    id_desi, d2d, _ = opt_coord[with_match].match_to_catalog_sky(desi_coord)
    match_gal = (d2d < spec_sep*u.arcsec)

    desi_z = np.ones(len(merged_all)) * np.nan
    desi_zconf = np.ones(len(merged_all)) * np.nan
    desi_sourceid = np.ones(len(merged_all), dtype='int') * -1
    
    desi_z[with_match[match_gal]] = desi['Z'][id_desi[match_gal]]
    desi_zconf[with_match[match_gal]] = desi['ZERR'][id_desi[match_gal]]
    desi_sourceid[with_match[match_gal]] = desi['TARGETID'][id_desi[match_gal]]

    merged_all.add_column(MaskedColumn(name='z_desi', data=desi_z,mask=(desi_z == np.nan), dtype=desi['Z'].dtype))
    merged_all.add_column(MaskedColumn(name='z_desi_err', data=desi_zconf,mask=(desi_zconf == np.nan), dtype=desi['ZERR'].dtype))
    merged_all.add_column(MaskedColumn(name='desi_sourceid', data=desi_sourceid,mask=(desi_sourceid == -1), dtype=desi['TARGETID'].dtype))
                                                                 
    """
    X-ray Matching
    """

    print('X-ray matching')

    rxs = Table.read(xray_2rxs_path)
    rxs_coord = SkyCoord(rxs['ALLW_RA'].data, rxs['ALLW_DEC'].data, unit='deg')

    id_rxs, d2d, _ = opt_coord[with_match].match_to_catalog_sky(rxs_coord)
    match_rxs = (d2d < xray_sep*u.arcsec)
    rxs_id = np.array(['']*len(opt_coord), dtype=rxs['2RXS_ID'].dtype)
    rxs_id[with_match[match_rxs]] = rxs['2RXS_ID'][id_rxs[match_rxs]]

    xmmsl2 = Table.read(xray_xmmsl2_path)
    xmmsl2_coord = SkyCoord(xmmsl2['ALLW_RA'].data, xmmsl2['ALLW_DEC'].data, unit='deg')

    id_xmm, d2d, _ = opt_coord[with_match].match_to_catalog_sky(xmmsl2_coord)
    match_xmm = (d2d < xray_sep*u.arcsec)

    xmm_id = np.array(['']*len(opt_coord), dtype=xmmsl2['XMMSL2_ID'].dtype)
    xmm_id[with_match[match_xmm]] = xmmsl2['XMMSL2_ID'][id_xmm[match_xmm]]

    merged_all.add_column(MaskedColumn(name='2RXS_ID', data=rxs_id,
                                       mask=(rxs_id=='')))
    merged_all.add_column(MaskedColumn(name='XMMSL2_ID', data=xmm_id,
                                       mask=(xmm_id=='')))


    t=merged_all

    t.write('temp.fits',overwrite=True)

    print('Remove -99s: ',end='')
    dblcols=[n for (n,ty) in t.dtype.descr if ('f8' in ty or 'f4' in ty)]
    for c in dblcols:
        print(c,end=' ')
        sys.stdout.flush()
        t[c]=np.where(t[c]==-99,np.nan,t[c])

    print('\nRemove 1e20s: ',end='')
    dblcols=[n for (n,ty) in t.dtype.descr if ('f8' in ty or 'f4' in ty)]
    for c in dblcols:
        print(c,end=' ')
        sys.stdout.flush()
        t[c]=np.where(t[c]==1e20,np.nan,t[c])

    print('\nRemove whitespace padding: ',end='')
    sys.stdout.flush()
    stringcols=[n for (n,ty) in t.dtype.descr if 'S' in ty]
    for c in stringcols:
        print(c,end=' ')
        sys.stdout.flush()
        t[c]=[str(s).rstrip() for s in t[c]]

    print('\nWrite to disk:')

    t.write(infile.replace('.fits',f'-withz-{hemisphere}.fits'), overwrite=True)

    final_tables.append(t)

if len(final_tables)>1:
    print('Stack and write joined table:')
    vstack(final_tables).write(infile.replace('.fits','-withz-joined.fits'), overwrite=True)

print('Done!')
