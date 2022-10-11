'''
add_gaus_info.py
add some info on the gaussians to the lofar source table
- it's a bit slow to do the separation calculation for all M sources, so do it once and save it
'''
import sys
import numpy as np
from astropy.table import Table, join, Column
import astropy.coordinates as ac


if __name__=='__main__':

    if len(sys.argv) == 1:
        print("Usage is : python add_gaus_info.py field_code ")
        print('E.g.: python add_gaus_info.py 0 ')
        sys.exit(1)

    h = str(sys.argv[1])
    if 'h' not in h:
        h+='h'
    if h not in  ['0h','13h','n0h','n13h','s0h','s13h']:
        print('unknown field code (should be 0h or 13h)',h)
        sys.exit(1)

    ### Required INPUTS
    # lofar source catalogue, gaussian catalogue and ML catalogues for each

    lLR_thresh = 0.404            # LR threshold
    
    redo = True

    path = '/Users/w.williams/projects/lofar_surveys/DR2/'
    lofargcat_file = path+'lr/LoTSS_DR2_{version}.gaus_{h}.lr-full.fits'.format(version=version,h=h)
    lofarcat_file = path+'LoTSS_DR2_{version}.srl_{h}.lr-full.presort.hdf5'.format(version=version,h=h)

    gaus_cols = ['Ng', 'G_max_sep', 'G_LR_max', 'Ng_LR_good','Ng_LR_good_unique','N_G_LR_matchsource','Flag_G_LR_problem']

    # Gaus catalogue
    lofargcat = Table.read(lofargcat_file)
    # only relevant gaussians are in M or C sources
    lofargcat = lofargcat[lofargcat['S_Code'] != 'S']

    # Source catalogue
    lofarcat = Table.read(lofarcat_file)

    if redo:
        for gaus_col in gaus_cols:
            if gaus_col in lofarcat.colnames:
                lofarcat.remove_column(gaus_col)

    lofarcat.add_column(Column(np.ones(len(lofarcat),dtype=int), 'Ng'))
    lofarcat.add_column(Column(np.ones(len(lofarcat),dtype=float), 'G_max_sep'))

    # gaus condensed LR info for msources
    lofarcat.add_column(Column(-1*np.ones(len(lofarcat),dtype=float), 'G_LR_max'))
    lofarcat.add_column(Column(-1*np.ones(len(lofarcat),dtype=int), 'Ng_LR_good'))
    lofarcat.add_column(Column(-1*np.ones(len(lofarcat),dtype=int), 'Ng_LR_good_unique'))
    lofarcat.add_column(Column(-1*np.ones(len(lofarcat),dtype=int), 'N_G_LR_matchsource'))
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=bool), 'Flag_G_LR_problem'))
    
    # Philips cases
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=int), 'G_LR_case1'))
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=int), 'G_LR_case2'))
    lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=int), 'G_LR_case3'))
    
    
    
    ## PS ML - matches for sources and gaussians
    #lofarcat = Table.read(lofarcat_file)
    #lofargcat = Table.read(lofargcat_file)
    
    # update to lr threshold but not yet catalogue
    fixlr = (lofarcat['lr'] < lLR_thresh)
    lofarcat['lr'][fixlr] = np.nan
    lofarcat['ra'][fixlr] = np.nan
    lofarcat['dec'][fixlr] = np.nan
    
    
    fixlr = (lofargcat['lr'] < lLR_thresh)
    lofargcat['lr'][fixlr] = np.nan
    lofargcat['ra'][fixlr] = np.nan
    lofargcat['dec'][fixlr] = np.nan
    
    ## renaming/fixing ML information
    lrcol = lofarcat['lr']
    lofarcat.add_column(Column(lrcol, 'cLR'))
    lrcol[np.isnan(lrcol)] = 0
    lofarcat.add_column(Column(lrcol, 'LR'))
    lrcol = np.zeros(len(lofarcat),dtype='S19')
    lrcol = lofarcat['UNWISE_OBJID']
    lofarcat.add_column(Column(lrcol, 'LR_name_wise'))
    lrcol = np.zeros(len(lofarcat),dtype='S19')
    lrcol = lofarcat['UID_L']
    lofarcat.add_column(Column(lrcol, 'LR_name_l'))
    lrcol = np.zeros(len(lofarcat),dtype=float)
    lrcol = lofarcat['ra']
    lofarcat.add_column(Column(lrcol, 'LR_ra'))
    lrcol = np.zeros(len(lofarcat),dtype=float)
    lrcol = lofarcat['dec']
    lofarcat.add_column(Column(lrcol, 'LR_dec'))


    lrgcol = lofargcat['lr']
    lofargcat.add_column(Column(lrgcol, 'LR'))
    lrgcol = np.zeros(len(lofargcat),dtype=float)
    lrgcol = lofargcat['ra']
    lofargcat.add_column(Column(lrgcol, 'LR_ra'))
    lrgcol = np.zeros(len(lofargcat),dtype=float)
    lrgcol = lofargcat['dec']
    lofargcat.add_column(Column(lrgcol, 'LR_dec'))
    lrgcol = np.zeros(len(lofargcat),dtype='S19')
    lrgcol = lofargcat['UNWISE_OBJID']
    lofargcat.add_column(Column(lrgcol, 'LR_name_wise'))
    lrgcol = np.zeros(len(lofargcat),dtype='S19')
    lrgcol = lofargcat['UID_L']
    lofargcat.add_column(Column(lrgcol, 'LR_name_l'))

    
    ## for single gaussians this is the same: we overwrite the M/C sources next
    #lofarcat.add_column(Column(lofarcat['lr_pc_7th'], 'gLR'))
    lrcol = lofarcat['lr']
    lofarcat.add_column(Column(lrcol, 'cgLR'))
    lrcol[np.isnan(lrcol)] = 0
    lofarcat.add_column(Column(lrcol, 'gLR'))
    lrcol = np.zeros(len(lofarcat),dtype='S19')
    lrcol = lofarcat['UNWISE_OBJID']
    lofarcat.add_column(Column(lrcol, 'gLR_name_wise'))
    lrcol = np.zeros(len(lofarcat),dtype='S19')
    lrcol = lofarcat['UID_L']
    lofarcat.add_column(Column(lrcol, 'gLR_name_l'))
    lrcol = np.zeros(len(lofarcat),dtype=float)
    lrcol = lofarcat['ra']
    lofarcat.add_column(Column(lrcol, 'gLR_ra'))
    lrcol = np.zeros(len(lofarcat),dtype=float)
    lrcol = lofarcat['dec']
    lofarcat.add_column(Column(lrcol, 'gLR_dec'))

    

    m_S = lofarcat['S_Code'] =='S'
    minds = np.where(~m_S)[0]
    print('Counting gaussians for {n} sources'.format(n=len(minds)))
    for i,sid in zip(minds, lofarcat['Source_Name'][~m_S]):
        ig = np.where(lofargcat['Source_Name']==sid)[0]
        Ng = len(ig)
        lofarcat['Ng'][i]= Ng
        if Ng == 0: ## this should not happen! but apparently does?
            print(sid, Ng)
            lofarcat['G_max_sep'][i] = np.nan
        else:
        
        
            gcoords = ac.SkyCoord(lofargcat['RA'][ig], lofargcat['DEC'][ig])
            _, sep, _ = gcoords.match_to_catalog_sky(gcoords, nthneighbor=Ng)
            lofarcat['G_max_sep'][i] = np.max(sep.to('arcsec').value)
        
        
        
    print('calculating LR stuff for m sources - this can be slow')
    m_S = lofarcat['S_Code'] =='S'
    #minds = np.where(~m_S)[0]
    minds = np.where(lofarcat['Ng']>1)[0]
    for i,sid in zip(minds, lofarcat['Source_Name'][~m_S]):
        ig = np.where(lofargcat['Source_Name']==sid)[0]
        lofarcat['G_LR_max'][i]= np.nanmax(lofargcat['lr'][ig])
        if np.any(np.isfinite(lofargcat['lr'][ig])):
            
            igi = np.nanargmax(lofargcat['lr'][ig])
            #for now, record the best gaussian
            #if lofarcat['G_LR_max'][i] > lofarcat['lr'][i]:
            lofarcat['gLR'][i] = lofarcat['G_LR_max'][i]
            lofarcat['gLR_name_l'][i] = lofargcat['LR_name_l'][ig[igi]]
            lofarcat['gLR_name_wise'][i] = lofargcat['LR_name_wise'][ig[igi]]
            lofarcat['gLR_ra'][i] = lofargcat['ra'][ig[igi]]
            lofarcat['gLR_dec'][i] = lofargcat['dec'][ig[igi]]
            #pass
        # how many unique acceptable matches are there for the gaussian components
        matches_ra = np.unique(lofargcat['ra'][ig][lofargcat['lr'][ig] >= lLR_thresh])
        matches_src_ra = (lofargcat['ra'][ig] == lofarcat['ra'][i]) & (lofargcat['dec'][ig] == lofarcat['dec'][i])
        lofarcat['N_G_LR_matchsource'][i] =  1*np.sum(matches_src_ra)
        n_matches_ra = len(matches_ra)
        if n_matches_ra > 1:
            lofarcat['Flag_G_LR_problem'][i] = True
        # any different to source match
        if np.sum(matches_ra != lofarcat['ra'][i]):
            lofarcat['Flag_G_LR_problem'][i] = True
        lofarcat['Ng_LR_good'][i]= np.nansum(lofargcat['lr'][ig] >= lLR_thresh)
        lofarcat['Ng_LR_good_unique'][i]= n_matches_ra
        
        
        # special case 1:
        # source with lr, 1 gaus with lr diff to source
        
        if (lofarcat['lr'][i] >= lLR_thresh) and (lofarcat['Ng_LR_good'][i] == 1) and (lofarcat['N_G_LR_matchsource'][i] == 0):
            #compare LRsource and LRgaus
            
            igi = ig[(lofargcat['lr'][ig] )>= lLR_thresh]
            if ((lofarcat['lr'][i] > 10 ) & (lofargcat['lr'][igi] < 10 ) & (lofarcat['lr'][i] > 10 * lofargcat['lr'][igi]))[0]:
                lofarcat['G_LR_case1'][i] = 1 # accept ML source
                #print '1'
            elif ((lofarcat['lr'][i] < 10 ) & (lofargcat['lr'][igi] > 10 ) & (lofargcat['lr'][igi] > 10 * lofarcat['lr'][i]))[0]:
                lofarcat['G_LR_case1'][i] = 2 # accept ML gaus
                #print '2'
            else:
                lofarcat['G_LR_case1'][i] = 3 # lgz
                #print '3'
            
        # special case 2:
        # source with lr, 2 gaus with lr, 1 diff to source
        
        elif (lofarcat['lr'][i] >= lLR_thresh) and (lofarcat['Ng_LR_good'][i] == 2) and (lofarcat['N_G_LR_matchsource'][i] == 1):
            # compare LRsource and LRgaus's
            igi = ig[(lofargcat['lr'][ig] )>= lLR_thresh]
            igs = igi[lofargcat['ra'][igi] == lofarcat['ra'][i]]
            ign = igi[lofargcat['ra'][igi] != lofarcat['ra'][i]]
            
            
            ## nb using max here not min as written in flowchart
            if ((lofarcat['lr'][i] > 100 ) and (lofargcat['lr'][igs] > 5*lofargcat['lr'][ign]) ):
                lofarcat['G_LR_case2'][i] = 1  # accept ML source
            elif (( lofarcat['lr'][i] < 30 ) and (lofargcat['lr'][igs] > np.max([10, lofarcat['lr'][i]]) ) and (lofargcat['lr'][ign] > np.max([10, lofarcat['lr'][i]]) ) ):
                lofarcat['G_LR_case2'][i] = 2  # deblend into 2 gaus with both ML
            else:
                lofarcat['G_LR_case2'][i] = 3  # lgz
        
        # special case 3:
        elif (lofarcat['lr'][i] < lLR_thresh) and (lofarcat['Ng_LR_good'][i] == 1):
            igi = ig[(lofargcat['lr'][ig] )>= lLR_thresh]
            if (lofargcat['lr'][igi] >  10*lLR_thresh) and (lofargcat['Maj'][igi] < 10):
                lofarcat['G_LR_case3'][i] = 1  # accept match
            else:
                lofarcat['G_LR_case3'][i] = 2  # lgz
    print('done')
                
        
        
    lofarcat['G_LR_max'][m_S] = lofarcat['lr'][m_S]
    lofarcat['Ng_LR_good'][m_S] = 1*(lofarcat['lr'][m_S] >= lLR_thresh)

    # some flags for mult_gaus sources:
    # source has good LR match, and no gaus
    # multiple matches to different sources
    # source has no good LR match, but one gaus does
        
        
        
    lofarcat.write(lofarcat_file, overwrite=True, serialize_meta=True)
        
