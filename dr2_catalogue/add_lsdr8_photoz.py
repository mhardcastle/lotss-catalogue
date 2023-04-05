# K Duncan - Jun 2021
# Cross-matching LoTSS DR2 Catalogues with LS DR8 Photo-z
#
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
dir=os.getcwd()
field=os.path.basename(dir)
print('field is',field)

g=sorted(glob.glob('final-v???.fits'))

infile=g[-1]
print('Processing',infile)

lofar_cat_path = infile
version_number = '0.3' # photoz merge version counter

spec_sep = 1.5 # arcsec - matching radius for spec-zs
xray_sep = 3 # arcsec - matching radius for X-ray cross-IDs

verbose = True

# work out what hemisphere(s) we're using

t = Table.read(lofar_cat_path)
if field=='Fall':
    critical_dec=36 # force use of South for all
else:
    critical_dec=32.275
    

hemispheres=[]

print('DEC range is',np.nanmin(t['ID_DEC']),np.nanmax(t['ID_DEC']))

if np.nanmax(t['ID_DEC'])>=critical_dec:
    hemispheres.append('north')
if np.nanmin(t['ID_DEC'])<critical_dec:
    hemispheres.append('south')

del(t)

final_tables=[]
for hemisphere in hemispheres:
    print('======= Doing hemisphere %s =======' % hemisphere)
    ancillary_data = '/beegfs/lofar/duncan/ancillary_data'

    photoz_cat_path = '/beegfs/lofar/duncan/photoz/{}_merged'.format(hemisphere)

    lofar_cat=Table.read(lofar_cat_path)
    filter_dec=np.where(np.isnan(lofar_cat['ID_DEC']),lofar_cat['DEC'],lofar_cat['ID_DEC'])
    
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
    xray_2rxs_path = '{0}/{1}'.format(ancillary_data, '2RXS_AllWISE_catalog_paper_2017May26.fits')
    xray_xmmsl2_path = '{0}/{1}'.format(ancillary_data, 'XMMSL2_AllWISE_catalog_paper_2017JUN09.fits')
    bricks_path = '{0}/{1}'.format(ancillary_data, 'survey-bricks.fits')

    """
    Pre-process LOFAR Catalogue to produce consistent Legacy Source IDs
    """

    # create dict to map brickid to brickname
    bricks = Table.read(bricks_path)
    brick_dict = dict(zip(bricks['BRICKNAME'], bricks['BRICKID']))
    brick_dict[''] = -1

    # Extract brick/object properties from UID
    ls_match = np.logical_and(lofar_cat['UID_L'] != '', lofar_cat['UID_L'] != 'N/A')
    brickname = np.array([(id.split('_')[0]) for id in lofar_cat['UID_L'][ls_match]])
    brickid = np.array([brick_dict[name] for name in brickname])
    objid = np.array([int(id.split('_')[1]) for id in lofar_cat['UID_L'][ls_match]])

    # Reformat into photoz-consistent Legacy ID
    newid = np.ones(len(lofar_cat), dtype='int')*-1
    newid[ls_match] = np.array(release*1e12 + brickid*1e6 + objid, dtype='int')

    lofar_cat['Legacy_ID'] = newid


    """
    Merge LR Matched sources through Legacy ID
    """

    with_match = (lofar_cat['ID_DEC'] > 0)

    #lofar_hpx_radio = healpy.ang2pix(2**3, lofar_cat['RA'], lofar_cat['DEC'],
    #                                 lonlat=True)

    # fake up co-ordinate set in order to pass non-null values to healpy
    # and allow every source to have a HPX -- only the ones with_match will be used

    ra=np.where(with_match,lofar_cat['ID_RA'],lofar_cat['RA'])
    dec=np.where(with_match,lofar_cat['ID_DEC'],lofar_cat['DEC'])

    lofar_hpx_opt = healpy.ang2pix(2**3, ra, dec, lonlat=True)

    lofar_cat['HPX']=lofar_hpx_opt

    lofar_hpx_unique = np.unique(lofar_hpx_opt)

    ls_joined_all = []

    # Loop through Healpix chunks to join photo-z information
    for hpx in lofar_hpx_unique[:]:
        if verbose:
            print(hpx)
        lofar_subset = (lofar_hpx_opt == hpx)
        hpx_opt = healpy.ang2pix(2**3, lofar_cat['ID_RA'][with_match*lofar_subset],
                                       lofar_cat['ID_DEC'][with_match*lofar_subset],
                                       lonlat=True)

        opt_hpx_unq = np.unique(lofar_hpx_opt)
        hpx_list = np.append(hpx, healpy.get_all_neighbours(2**3, hpx))
        hpx_list = hpx_list[np.in1d(hpx_list, opt_hpx_unq)] # Crop surplus neighbours
        #print(opt_hpx_unq)

        pzpath = '{0}/hpx_{1:03}_merged.fits'
        filenames=[pzpath.format(photoz_cat_path, o) for o in hpx_list]
        if verbose:
            print(filenames)
        tables=[Table.read(f) for f in filenames if os.path.isfile(f)]
        if len(tables):
            print(len(tables))
            try:
                photoz = vstack(tables)
            except TypeError:
                print(tables)
                raise

            photoz['id'].name = 'Legacy_ID'

            ls_joined = join(lofar_cat[lofar_subset], photoz,
                             keys='Legacy_ID', join_type='left')

            if ls_joined is not None:
                ls_joined_all.append(ls_joined)


    if not ls_joined_all:
        continue
    try:
        merged_all = vstack(ls_joined_all)
    except TypeError:
        print('Failed on',ls_joined_all)
        raise
    merged_all.sort('Source_Name')

    with_match = np.where(merged_all['ID_DEC'] > 0)[0]
    opt_coord = SkyCoord(merged_all['ID_RA'], merged_all['ID_DEC'], unit='deg')

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

    print('Blank columns where ID does not match: ',end='')
    filt=(t['Legacy_ID']==-1)
    filt|=t['release'].mask
    cols=['ra','dec','pstar','zphot','zphot_err']
    for c in t.colnames:
        if 'fracflux' in c: cols.append(c)
        if 'lupt' in c: cols.append(c)
        if 'var.' in c: cols.append(c)
    for c in cols:
        print(c,end=' ')
        sys.stdout.flush()
        t[c]=np.where(filt,np.nan,t[c])

    #t['Legacy_ID']=np.where(filt,np.nan,t['Legacy_ID'])

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

    t.write('{0}_photoz_v{1}_{2}.fits'.format(os.path.splitext(os.path.split(lofar_cat_path)[1])[0],version_number,hemisphere), overwrite=True)

    final_tables.append(t)

if len(final_tables)>1:
    print('Stack and write joined table:')
    vstack(final_tables).write('{0}_photoz_v{1}_{2}.fits'.format(os.path.splitext(os.path.split(lofar_cat_path)[1])[0],version_number,'joined'), overwrite=True)

print('Done!')
