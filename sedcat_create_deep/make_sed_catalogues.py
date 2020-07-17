import sys
import os
from astropy.table import Table, Column, vstack
import astropy.units as au
import numpy as np
import utils.plot_util as pp

#field = 'Bootes'
#field = 'LH'
#field = 'EN1'
field = sys.argv[1]

if field=='all':
    fields = ['Bootes','EN1','LH']
else:
    fields = [field]


version = 'v0.7_newtest'
version = 'v0.8'
version = 'v1.0'

zpapply = 'atlas'


for field in fields:
    path = '/beegfs/lofar/wwilliams/lofar_surveys/deep/science_ready_catalogs/filter_information/'
    filterpathFIR = 'FIR_filters_filters/'
    if field == 'EN1':
        finfo = path+'EN1_filters.res.info'
        filterpath = 'EN1_filters_filters/'
        zpfile = path+'../zeropoint_offsets/en1_merged_zeropoint_offsets.txt'
        outfilterpathherschel = 'filters/Herschel/'
        ftrans = path+'EN1.filter.translate'
        datfile = 'agnfitter_cols_en1_wl.csv'
        
        fxmatchin = '/beegfs/lofar/deepfields/lgz/en1/final-v1.0.fits'

        fchangev0p8id = '/beegfs/lofar/deepfields/changes_v0.7_v0.8/EN1_v0.6_v0.7_changedIDs_v0.8.fits'
        fchangev0p8z = '/beegfs/lofar/deepfields/changes_v0.7_v0.8/EN1_v0.7_v0.8_changedZBEST.fits'
        fchangev1p0fir = '/beegfs/lofar/deepfields/changes_v0.7_v0.8/EN1_changedFIRfluxes_v0.8_v1.0.fits'
        
              
        # outlier info from:
        ftoutlier = '/beegfs/lofar/deepfields/science_ready_catalogs/EN1_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits'
        
        
    elif field == 'Bootes':
        finfo = path+'filter.bootes_mbrown_2014a.res.info'
        filterpath = 'filter.bootes_mbrown_2014a_filters/'
        zpfile = path+'../zeropoint_offsets/bootes_merged_zeropoint_offsets.txt'
        outfilterpathherschel = 'filters/Herschel/'
        ftrans = path+'brown.zphot.2014.translate'
        datfile = 'agnfitter_cols_bootes_wl.csv'
        
        fxmatchin = '/beegfs/lofar/deepfields/lgz/bootes/final-v1.0.fits'
        
        
        fchangev0p8id = '/beegfs/lofar/deepfields/changes_v0.7_v0.8/Bootes_v0.6_v0.7_changedIDs_v0.8.fits'
        fchangev1p0fir = '/beegfs/lofar/deepfields/changes_v0.7_v0.8/Bootes_changedFIRfluxes_v0.7_v1.0.fits'
        fchangev1p0fd = '/beegfs/lofar/deepfields/changes_v0.7_v0.8/Bootes_changedFLAGDEEP_v0.7_v1.0.fits'

                
        # outlier info from:
        ftoutlier = '/beegfs/lofar/deepfields/science_ready_catalogs/Bootes_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits'
        
        
    elif field == 'LH':
        finfo = path+'Lockman-SWIRE_filters.res.info'
        filterpath = 'Lockman-SWIRE_filters_filters/'
        outfilterpathherschel = 'filters/Herschel/'
        zpfile = path+'../zeropoint_offsets/lh_merged_zeropoint_offsets.txt'
        ftrans = path+'LH.filter.translate'
        datfile = 'agnfitter_cols_lockman_wl.csv'
        
        fxmatchin = '/beegfs/lofar/deepfields/lgz/lockman/final-v1.0.fits'  # same as on surveys webpage???
        
        
        fchangev0p8id = '/beegfs/lofar/deepfields/changes_v0.7_v0.8/Lockman_v0.6_v0.7_changedIDs_v0.8.fits'
        fchangev0p8z = '/beegfs/lofar/deepfields/changes_v0.7_v0.8/Lockman_v0.7_v0.8_changedZBEST.fits'
        fchangev1p0fir = '/beegfs/lofar/deepfields/changes_v0.7_v0.8/Lockman_changedFIRfluxes_v0.8_v1.0.fits'
        
        # outlier info from:
        ftoutlier = '/beegfs/lofar/deepfields/science_ready_catalogs/LH_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits'
        
        
    else:
        print('not implemented')
        sys.exit(1)
        
        
    ### merge new fir 
    
    txmatchin = Table.read(fxmatchin)
    
    
    #### should be merged in already since v0.7 ###
    #colsnew = tnew.colnames
    #cols = txmatchin.colnames

    #copycols = []
    #for cc in colsnew:
        #if cc in cols:
            #copycols.append(cc)

    ## for the sources in tnew - replace the columns
    #minds = []
    #for ti in tnew:
        #mind = np.where(txmatchin['Source_Name'] == ti['Source_Name'])[0][0]
        #minds.append(mind)
        #print(ti['Source_Name'],mind)
    #for cc in copycols:
        #txmatchin[cc][minds] = tnew[cc]


    ### add a flag for these sources...
    #newflag = np.zeros(len(t), dtype=bool)
    #newflag[minds] = 1
    #txmatchin.add_column(Column(name='newXIDp', data=newflag))


    ## add SERVS moc info for field EN1 only
    #if field == 'EN1':
        #print('get moc flags')
        ## for the sources in tmoc - add the column
        #minds = []
        #mindst = []
        #for ti in txmatchin:
            #if (ti['ID'] in tmoc['ID']):
                #mind = np.where(ti['ID'] == tmoc['ID'])[0][0]
                #minds.append(mind)
                #mindst.append(True)
                ##print ti['Source_Name'],mind
            #else:
                #mindst.append(False)
            
        #txmatchin.add_column(Column(name='FLAG_SERVS', data=np.zeros(len(txmatchin),dtype=bool)))
        #txmatchin['FLAG_SERVS'] = tmoc['FLAG_SERVS'][minds]

    def get_matchind(t1, t2):
        inds = []
        for s in t2['Source_Name']:
            inds.append(np.where(t1['Source_Name'] == s)[0][0])
        return inds


    ## add flags if the source appears in the changes files (recording changes to v0.7 and to v0.8
    print('get version flags')
    
    
    txmatchin.add_column(Column(name='CHANGE_FLAG_ID', data=np.zeros(len(txmatchin),dtype=bool)))
    tchangev0p8id = Table.read(fchangev0p8id)
    inds = get_matchind(txmatchin, tchangev0p8id)
    txmatchin['CHANGE_FLAG_ID'][inds] = True
    
    
    txmatchin.add_column(Column(name='CHANGE_FLAG_ZBEST', data=np.zeros(len(txmatchin),dtype=bool)))
    txmatchin.add_column(Column(name='CHANGE_FLAG_DEEP', data=np.zeros(len(txmatchin),dtype=bool)))
    if field != 'Bootes':

        tchangev0p8z = Table.read(fchangev0p8z)
        inds = get_matchind(txmatchin, tchangev0p8z)
        txmatchin['CHANGE_FLAG_ZBEST'][inds] = True
    else:
        
        tchangev1p0fd = Table.read(fchangev1p0fd)
        inds = get_matchind(txmatchin, tchangev1p0fd)
        txmatchin['CHANGE_FLAG_DEEP'][inds] = True
        
    txmatchin.add_column(Column(name='CHANGE_FLAG_FIR', data=np.zeros(len(txmatchin),dtype=bool)))
    tchangev1p0fir = Table.read(fchangev1p0fir)
    inds = get_matchind(txmatchin, tchangev1p0fir)
    txmatchin['CHANGE_FLAG_FIR'][inds] = True
        



    ftfir = 'temp_mergenew.fits'
    # save so we can stilts match with outliers
    txmatchin.write(ftfir,overwrite=True)

    #### done with new XID+ merge



    #### get outlier flags
    print('getting outlier flags')
    ftmatch = 'temp_matchout.fits'
    
    os.system('stilts tmatch2 in1={in1} in2={in2} out={out} values1=ID values2=id join=all1 find=best suffix2=outlier matcher=exact'.format(in1=ftfir, in2=ftoutlier, out=ftmatch))
    


    cat = Table.read(ftfir)
    toutlier = Table.read(ftmatch)

    outcat = Table.read(ftfir)
    outcat.keep_columns(['ID'])
        
    outcat.write(field+'_test_{version}.fits'.format(version=version), overwrite=True)

    '''
    ### Bootes
    ----------

    FLAG_CLEAN == 1
    FLAG_DEEP != 0

    I_fluxerr > 0
    ch2_fluxerr > 0

    ### ELAIS-N1
    ------------

    FLAG_CLEAN == 1

    i_fluxerr > 0
    ch2_swire_fluxerr > 0

    or

    (Note, the FLAG_OVERLAP == 7 being applied for LR matching will supercede these cuts)

    ### Lockman-Hole
    ----------------

    FLAG_CLEAN == 1

    r_fluxerr > 0
    ch2_swire_fluxerr > 0
    '''
        
    #toutlier.add_column(Column(np.zeros(len(toutlier),dtype=bool),'FLAG_GOOD'))
    
    cat.add_column(Column(np.zeros(len(cat),dtype=bool),'FLAG_GOOD'))
    if field == 'EN1':
        cat['FLAG_GOOD'][(np.isfinite(toutlier['ID_1'])) & (toutlier['FLAG_CLEANoutlier'] == 1) & (toutlier['i_fluxerr']>0) & (toutlier['ch2_swire_fluxerr']>0)] = True
    elif field == 'Bootes':
        cat['FLAG_GOOD'][(np.isfinite(toutlier['ID_1'])) &(toutlier['FLAG_CLEANoutlier'] == 1) & (toutlier['FLAG_DEEPoutlier'] != 0) & (toutlier['I_fluxerr']>0) & (toutlier['ch2_fluxerr']>0)] = True
    elif field == 'LH':
        cat['FLAG_GOOD'][(np.isfinite(toutlier['ID_1'])) &(toutlier['FLAG_CLEANoutlier'] == 1) & (toutlier['r_fluxerr']>0) & (toutlier['ch2_swire_fluxerr']>0)] = True
    


    #sys.exit()

    for col in cat.colnames:
        if 'flux_corr' in col:
            col1 = col.replace('_corr','')
            if col1 not in toutlier.colnames:  # handle the new FUV/NUV that aren't in Ken's files
                print(col1)
                continue
            Table([cat[col], toutlier[col1]])
            print( col1, col, np.sum(toutlier[col1] == -90),np.sum(np.isnan(cat[col])),)
            cat[col][toutlier[col1] == -90] = np.nan
            print(np.sum(np.isnan(cat[col])))

    #### done with outliers



    outfilterpath = 'filters/'+field+'/'

    # filter information is stored in this file
    filterdat = Table.read('filters/'+datfile)
    
    ## add GALEX FUV and NUV filters
    if field != 'Bootes':
        datfile_bootes = 'agnfitter_cols_bootes_wl.csv'
        filterdat_bootes = Table.read('filters/'+datfile_bootes)
        
        filterdat = vstack((filterdat_bootes[filterdat_bootes['id'] =='FUV'], filterdat))
        filterdat = vstack((filterdat_bootes[filterdat_bootes['id'] =='NUV'], filterdat))
        
    ## flag to include or exclude FUV and NUV for all
    stripUV = True
    if stripUV:
        filterdat = filterdat[filterdat['id'] != 'FUV']
        filterdat = filterdat[filterdat['id'] != 'NUV']
    
    
    if not os.path.exists(outfilterpath):
        os.mkdir(outfilterpath)




    cat.add_column(Column(np.arange(len(cat)),'radioID'))
   
    
    # apply zeropoint offsets - multiply fluxes
    zpoffsets = Table.read(zpfile, format='ascii')
    for zpi in zpoffsets:
        zfilt = zpi['filter'] +'_corr' ## cat actually named filt_flux_corr
        if zfilt in cat.colnames:
            cat[zfilt] = cat[zfilt]*zpi[zpapply]
            
            print('scaling ',zfilt,'by', zpi[zpapply])
        else:
            print('no column:',zfilt)



    outcolnames = ['Source_Name', 'radioID', 'ID', 'z1_median', 'Z_BEST', 'XID+_rerun_mips', 'XID+_rerun_pacs', 'XID+_rerun_SPIRE', 'FLAG_GOOD']
    for cflag in ['CHANGE_FLAG_ID', 'CHANGE_FLAG_ZBEST', 'CHANGE_FLAG_FIR', 'CHANGE_FLAG_DEEP']:
        if cflag in cat.colnames:
            outcolnames.append(cflag)
    if 'FLAG_CLEAN' in cat.colnames:
        outcolnames.append('FLAG_CLEAN')
    if 'FLAG_DEEP' in cat.colnames:
        outcolnames.append('FLAG_DEEP')
    if 'FLAG_OVERLAP' in cat.colnames:
        outcolnames.append('FLAG_OVERLAP')
    #outcols = [cat[c] for c in outcolnames]
    fcat = cat.copy()
    cat.keep_columns(outcolnames)
        
        
    filterlist = []
    for fi in range(len(filterdat)):

        fluxcol = filterdat['influx_colname'][fi]
        fluxcolname = filterdat['outflux_colname'][fi]
        errcol = filterdat['influxerr_colname'][fi]
        errcolu = filterdat['influxerru_colname'][fi]
        errcoll = filterdat['influxerrl_colname'][fi]
        errcolname = filterdat['outfluxerr_colname'][fi]
        filt =  filterdat['id'][fi]
        
        filterlist.append(filterdat['id'][fi])
        fluxes = fcat[fluxcol]
        fluxes[fluxes==-99] = np.nan
        fluxes[fluxes==-90] = np.nan
        if fluxes.unit=='muJy':
            fluxes.unit='microJansky'
        ### mips is unitless but is in mJy (maybe)
        if fluxes.unit is None:
            
            if ('PACS' in fluxcolname) or ('SPIRE' in fluxcolname):
                fluxes = fluxes *au.mJy
                print('using mJy for',filt)
            #elif 'MIPS' in fluxcolname:
                #fluxes = fluxes *au.microJansky
                #print 'using microJansky for',filt
            else:
                #print 'not implemented', fluxcolname
                fluxes = fluxes *au.microJansky
                print('using muJy for',filt)
                #sys.exit()
        fluxes = fluxes.to(au.Jy)
        if errcol != '':
            errfluxes = fcat[errcol]
            errfluxes[errfluxes==-99] = np.nan
            errfluxes[errfluxes==-90] = np.nan
            
            if errfluxes.unit=='muJy':
                errfluxes.unit='microJansky'
            ### mips is unitless but is in mJy (maybe)
            if errfluxes.unit is None:
                
                if ('PACS' in fluxcolname) or ('SPIRE' in fluxcolname):
                    errfluxes = errfluxes *au.mJy
                    print('using mJy for',filt)
                #elif 'MIPS' in fluxcolname:
                    #errfluxes = errfluxes *au.microJansky
                    #print 'using microJansky for',filt
                else:
                    #print 'not implemented', fluxcolname
                    errfluxes = errfluxes *au.microJansky
                    print('using muJy for',filt)
                    #sys.exit()
            errfluxes = errfluxes.to(au.Jy)
        else:
            errfluxesl = fcat[errcolu]
            errfluxesu = fcat[errcoll]
            errfluxesl[errfluxesl==-99] = np.nan
            errfluxesu[errfluxesu==-99] = np.nan
            errfluxesl[errfluxesl==-90] = np.nan
            errfluxesu[errfluxesu==-90] = np.nan
            if errfluxesl.unit=='muJy':
                errfluxesl.unit='microJansky'
            if errfluxesl.unit is None:
                if ('PACS' in fluxcolname) or ('SPIRE' in fluxcolname):
                    errfluxesl = errfluxesl *au.mJy
                    print('using mJy for',filt)
                #elif 'MIPS' in fluxcolname:
                    #errfluxesl = errfluxesl *au.microJansky
                    #print 'using microJansky for',filt
                else:
                    #print 'not implemented', fluxcolname
                    errfluxesl = errfluxesl *au.microJansky
                    print('using muJansky for',filt)
            if errfluxesu.unit=='muJy':
                errfluxesu.unit='microJansky'
            if errfluxesu.unit is None:
                if ('PACS' in fluxcolname) or ('SPIRE' in fluxcolname):
                    errfluxesu = errfluxesu *au.mJy
                    print('using mJy for',filt)
                #elif 'MIPS' in fluxcolname:
                    #errfluxesu = errfluxesu *au.microJansky
                    #print 'using microJansky for',filt
                else:
                    #print 'not implemented', fluxcolname
                    errfluxesu = errfluxesu *au.microJansky
                    print('using microJansky for',filt)
            errfluxesl = errfluxesl.to(au.Jy)
            errfluxesu = errfluxesu.to(au.Jy)
            errfluxesl = np.abs(fluxes-errfluxesl)
            errfluxesu = np.abs(errfluxesu-fluxes)
            errfluxes = np.nanmax((errfluxesl, errfluxesu),axis=0) *au.Jy
            
        #### NB ###
        # add 10% flux in quadrature to flux errors
        errfluxes_raw = errfluxes.copy()
        errfluxes = np.sqrt(errfluxes**2. + (0.1*fluxes)**2. )
            
        # save columns
        cat.add_column(Column(fluxes, fluxcolname, unit=au.Jy))
        cat.add_column(Column(errfluxes, errcolname, unit=au.Jy))
        
        
        # keep flux errors without the 10% errors
        cat.add_column(Column(errfluxes_raw, errcolname+'_raw', unit=au.Jy))
        
    #outcat = Table(outcols)
    cat.write(field+'_sedfit_{version}.fits'.format(version=version), overwrite=True)
    
    
    ### No longer needed for new version
    ### add _wl cols for agnfitter
    for fi in range(len(filterdat)):
        #fcol = filterdat['fluxcol'][fi]
        #if fcol == '': continue

        #fluxcol = filterdat['fluxcol'][fi]
        lamcol = filterdat['outflux_colname'][fi].replace('_f','_wl')
        lam = filterdat['af_wavelength'][fi]
        cat.add_column(Column(lam*np.ones(len(cat)), lamcol, unit=au.angstrom))
        
    # remove the raw flux error columns from agnfitter
    outcols = cat.colnames
    outcols = [c for c in outcols if '_raw' not in c]
    cat.keep_columns(outcols)
    #outcat = Table(outcols)
    
    
    
    print('unique wavelengths', len(np.unique(filterdat['af_wavelength'])) )
    print('wavelengths', len(filterdat['af_wavelength']>0))
    if  len(np.unique(filterdat['af_wavelength']))  != len(filterdat['af_wavelength']):
        print('duplicate wavelengths')


    cat.write(field+'_agnfitter_{version}.fits'.format(version=version), overwrite=True)

    #Hfilternames = {1: 'PACS_100mu.txt', 2:  'PACS_160mu.txt',  3:'SPIRE_250mu.txt', 4:'SPIRE_350mu.txt', 5:'SPIRE_500mu.txt'}
    with open(field+'_'+version+'_agnfitter_filters.lst','w') as f:
        f.write('name|path|-|wavelength|\n')
        for fi in range(len(filterdat)):
            ffilter = filterdat['af_filtername'][fi]
            # update - new agnfitter...
            f.write(filterdat['outflux_colname'][fi].replace('_f','')+'|'+outfilterpath+ffilter+'.filter'+'|-|'+str(filterdat['af_wavelength'][fi])+'\n')
        


    # create new settings file
    settingsin = '/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/SETTINGS_AGNfitter_deep_template.py'
    settingsout = settingsin.replace('template',field+'_'+version)
    with open(settingsin,'r') as fin:
        with open(settingsout,'w') as fout:
            for line in fin.readlines():
                if '##FIELD##' in line:
                    line = line.replace('##FIELD##',field)
                if '##VERSION##' in line:
                    line = line.replace('##VERSION##',version)
                fout.write(line)


