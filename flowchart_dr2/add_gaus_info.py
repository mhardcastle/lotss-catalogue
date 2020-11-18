'''
add_gaus_info.py
add some info on the gaussians to the lofar source table
- it's a bit slow to do the separation calculation for all M sources, so do it once and save it
'''
import sys
from tqdm import tqdm
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

    
    ##LR threshold
    ##s0h -> 0.394
    ##s13h -> 0.328
    ##n13h -> 0.309  -- note, adopt this level for both N and S

    if h == '13h':
        lLR_thresh_n = 0.309            # LR threshold
        lLR_thresh_s = 0.328            # LR threshold
        LR_thresh_dec = 32.375
    elif h == '0h':
        lLR_thresh_n = 0.394            # LR threshold
        lLR_thresh_s = 0.394            # LR threshold
        LR_thresh_dec = -90
    else:
        print('LR threshold not implemented for field',h)
        sys.exit(1)
    
    redo = True

    path = '/data2/wwilliams/projects/lofar_surveys/LoTSS-DR2-Feb2020/'
    lofargcat_file = path+'lr/LoTSS_DR2_v100.gaus_{h}.lr-full.fits'.format(h=h)
    lofarcat_file = path+'LoTSS_DR2_v100.srl_{h}.lr-full.presort.hdf5'.format(h=h)

    gaus_cols = ['Ng', 'G_max_sep', 'G_LR_max', 'Ng_LR_good', 'Ng_LR_good_unique', 'N_G_LR_matchsource', 'Flag_G_LR_problem', 'G_LR_case1', 'G_LR_case2', 'G_LR_case3',                 'cLR', 'LR_threshold', 'LR_threshold10', 'LR', 'LR_name_wise', 'LR_name_l', 'LR_ra', 'LR_dec', 'gLR', 'gLR_name_wise','gLR_name_l','gLR_ra','gLR_dec']
    #gaus_cols = ['Ng', 'G_max_sep', 'G_LR_max', 'Ng_LR_good','Ng_LR_good_unique','N_G_LR_matchsource','Flag_G_LR_problem', 
                 #'G_LR_case1','G_LR_case2','G_LR_case3',
                 #'cLR','LR_threshold','LR_threshold10','LR','LR_name_wise','LR_name_l','LR_ra','LR_dec',
                 #'cgLR','gLR','gLR_name_wise','gLR_name_l','gLR_ra','gLR_dec']

    # Gaus catalogue
    lofargcat = Table.read(lofargcat_file)
    # only relevant gaussians are in M or C sources
    lofargcat = lofargcat[lofargcat['S_Code'] != 'S']

    # Source catalogue
    lofarcat = Table.read(lofarcat_file)
    Nlofarcat = len(lofarcat)
    Nlofargcat = len(lofargcat)

    if redo:
        for gaus_col in gaus_cols:
            if gaus_col in lofarcat.colnames:
                lofarcat.remove_column(gaus_col)

    
    lofarcat.add_column(Column(np.ones(Nlofarcat,dtype=int), 'Ng'))
    lofarcat.add_column(Column(np.ones(Nlofarcat,dtype=float), 'G_max_sep'))

    # gaus condensed LR info for msources
    lofarcat.add_column(Column(-1*np.ones(Nlofarcat,dtype=float), 'G_LR_max'))
    lofarcat.add_column(Column(-1*np.ones(Nlofarcat,dtype=int), 'Ng_LR_good'))
    lofarcat.add_column(Column(-1*np.ones(Nlofarcat,dtype=int), 'Ng_LR_good_unique'))
    lofarcat.add_column(Column(-1*np.ones(Nlofarcat,dtype=int), 'N_G_LR_matchsource'))
    lofarcat.add_column(Column(np.zeros(Nlofarcat,dtype=bool), 'Flag_G_LR_problem'))
    
    # Philips cases
    lofarcat.add_column(Column(np.zeros(Nlofarcat,dtype=int), 'G_LR_case1'))
    lofarcat.add_column(Column(np.zeros(Nlofarcat,dtype=int), 'G_LR_case2'))
    lofarcat.add_column(Column(np.zeros(Nlofarcat,dtype=int), 'G_LR_case3'))
    
    
    for cat in [lofarcat, lofargcat]:
        Ncat = len(cat)
        cat.add_column(Column(np.zeros(Ncat,dtype=bool), 'LR_threshold'))
        cat['LR_threshold'][(cat['DEC']>=LR_thresh_dec)&(cat['lr'] >=lLR_thresh_n)] = True
        cat['LR_threshold'][(cat['DEC']<LR_thresh_dec)&(cat['lr'] >=lLR_thresh_s)] = True
        cat.add_column(Column(np.zeros(Ncat,dtype=bool), 'LR_threshold10'))
        cat['LR_threshold10'][(cat['DEC']>=LR_thresh_dec)&(cat['lr'] >=10*lLR_thresh_n)] = True
        cat['LR_threshold10'][(cat['DEC']<LR_thresh_dec)&(cat['lr'] >=10*lLR_thresh_s)] = True
    
    
        mthresh = cat['LR_threshold']
        cat.add_column(Column(np.nan*np.ones(Ncat,dtype=bool), 'LR'))
        cat['LR'][mthresh] = cat['lr'][mthresh]
        cat.add_column(Column(np.zeros(Ncat,dtype='S19'), 'LR_name_wise'))
        cat['LR_name_wise'][mthresh] = cat['UNWISE_OBJID'][mthresh]   
        cat.add_column(Column(np.zeros(Ncat,dtype='S19'), 'LR_name_l'))
        cat['LR_name_l'][mthresh] = cat['UID_L'][mthresh]    
        cat.add_column(Column(np.nan*np.ones(Ncat,dtype=float), 'LR_ra'))
        cat['LR_ra'][mthresh] = cat['ra'][mthresh]    
        cat.add_column(Column(np.nan*np.ones(Ncat,dtype=float), 'LR_dec'))
        cat['LR_dec'][mthresh] = cat['dec'][mthresh]    
    

    
    ### for single gaussians this is the same: we overwrite the M/C sources next
    ##lofarcat.add_column(Column(lofarcat['lr_pc_7th'], 'gLR'))
    lrcol = lofarcat['LR']
    lofarcat.add_column(Column(lrcol, 'gLR'))
    lrcol = np.zeros(Nlofarcat,dtype='S19')
    lrcol = lofarcat['UNWISE_OBJID']
    lofarcat.add_column(Column(lrcol, 'gLR_name_wise'))
    lrcol = np.zeros(Nlofarcat,dtype='S19')
    lrcol = lofarcat['UID_L']
    lofarcat.add_column(Column(lrcol, 'gLR_name_l'))
    lrcol = np.zeros(Nlofarcat,dtype=float)
    lrcol = lofarcat['ra']
    lofarcat.add_column(Column(lrcol, 'gLR_ra'))
    lrcol = np.zeros(Nlofarcat,dtype=float)
    lrcol = lofarcat['dec']
    lofarcat.add_column(Column(lrcol, 'gLR_dec'))

    

    m_S = lofarcat['S_Code'] =='S'
    minds = np.where(~m_S)[0]
    print('Counting gaussians for {n} sources'.format(n=len(minds)))
    print('calculating LR stuff for m sources - this can be slow')
    for i,sid in tqdm(zip(minds, lofarcat['Source_Name'][~m_S])):
        ig = np.where(lofargcat['Source_Name']==sid)[0]
        Ng = len(ig)
        lofarcat['Ng'][i]= Ng
        if Ng == 0: ## this should not happen! but apparently does? it happens for S_Code=C
            print(sid, Ng)
            lofarcat['G_max_sep'][i] = 0.
        else:
        
        
            gcoords = ac.SkyCoord(lofargcat['RA'][ig], lofargcat['DEC'][ig])
            _, sep, _ = gcoords.match_to_catalog_sky(gcoords, nthneighbor=Ng)
            lofarcat['G_max_sep'][i] = np.max(sep.to('arcsec').value)
        
        
        
    #for i,sid in tqdm(zip(minds, lofarcat['Source_Name'][~m_S])):
        #ig = np.where(lofargcat['Source_Name']==sid)[0]
        lofarcat['G_LR_max'][i]= np.nanmax(lofargcat['LR'][ig])
        if np.any(np.isfinite(lofargcat['LR'][ig])):
            
            igi = np.nanargmax(lofargcat['LR'][ig])
            #for now, record the best gaussian
            #if lofarcat['G_LR_max'][i] > lofarcat['LR'][i]:
            lofarcat['gLR'][i] = lofarcat['G_LR_max'][i]
            lofarcat['gLR_name_l'][i] = lofargcat['LR_name_l'][ig[igi]]
            lofarcat['gLR_name_wise'][i] = lofargcat['LR_name_wise'][ig[igi]]
            lofarcat['gLR_ra'][i] = lofargcat['ra'][ig[igi]]
            lofarcat['gLR_dec'][i] = lofargcat['dec'][ig[igi]]
            #pass
        # how many unique acceptable matches are there for the gaussian components
        matches_ra = np.unique(lofargcat['ra'][ig][lofargcat['LR_threshold'][ig]])
        matches_src_ra = (lofargcat['ra'][ig] == lofarcat['ra'][i]) & (lofargcat['dec'][ig] == lofarcat['dec'][i])
        lofarcat['N_G_LR_matchsource'][i] =  1*np.sum(matches_src_ra)
        n_matches_ra = len(matches_ra)
        if n_matches_ra > 1:
            lofarcat['Flag_G_LR_problem'][i] = True
        # any different to source match
        if np.sum(matches_ra != lofarcat['ra'][i]):
            lofarcat['Flag_G_LR_problem'][i] = True
        lofarcat['Ng_LR_good'][i]= np.nansum(lofargcat['LR_threshold'][ig])
        lofarcat['Ng_LR_good_unique'][i]= n_matches_ra
        
        
        # special case 1:
        # source with lr, 1 gaus with lr diff to source
        
        if (lofarcat['LR_threshold'][i]) and (lofarcat['Ng_LR_good'][i] == 1) and (lofarcat['N_G_LR_matchsource'][i] == 0):
            #compare LRsource and LRgaus
            
            igi = ig[lofargcat['LR_threshold'][ig]]
            if ((lofarcat['LR'][i] > 10 ) & (lofargcat['LR'][igi] < 10 ) & (lofarcat['LR'][i] > 10 * lofargcat['LR'][igi]))[0]:
                lofarcat['G_LR_case1'][i] = 1 # accept ML source
                #print '1'
            elif ((lofarcat['LR'][i] < 10 ) & (lofargcat['LR'][igi] > 10 ) & (lofargcat['LR'][igi] > 10 * lofarcat['LR'][i]))[0]:
                lofarcat['G_LR_case1'][i] = 2 # accept ML gaus
                #print '2'
            else:
                lofarcat['G_LR_case1'][i] = 3 # lgz
                #print '3'
            
        # special case 2:
        # source with lr, 2 gaus with lr, 1 diff to source
        
        elif (lofarcat['LR_threshold'][i]) and (lofarcat['Ng_LR_good'][i] == 2) and (lofarcat['N_G_LR_matchsource'][i] == 1):
            # compare LRsource and LRgaus's
            igi = ig[lofargcat['LR_threshold'][ig]]
            igs = igi[lofargcat['ra'][igi] == lofarcat['ra'][i]]
            ign = igi[lofargcat['ra'][igi] != lofarcat['ra'][i]]
            
            
            ## nb using max here not min as written in flowchart
            if ((lofarcat['LR'][i] > 100 ) and (lofargcat['LR'][igs] > 5*lofargcat['LR'][ign]) ):
                lofarcat['G_LR_case2'][i] = 1  # accept ML source
            elif (( lofarcat['LR'][i] < 30 ) and (lofargcat['LR'][igs] > np.max([10, lofarcat['LR'][i]]) ) and (lofargcat['LR'][ign] > np.max([10, lofarcat['LR'][i]]) ) ):
                lofarcat['G_LR_case2'][i] = 2  # deblend into 2 gaus with both ML
            else:
                lofarcat['G_LR_case2'][i] = 3  # lgz
        
        # special case 3:
        elif (~lofarcat['LR_threshold'][i]) and (lofarcat['Ng_LR_good'][i] == 1):
            igi = ig[lofargcat['LR_threshold'][ig] ]
            if (lofargcat['LR_threshold10'][igi] ) and (lofargcat['Maj'][igi] < 10):
                lofarcat['G_LR_case3'][i] = 1  # accept match
            else:
                lofarcat['G_LR_case3'][i] = 2  # lgz
    print('done')
                
        
        
    lofarcat['G_LR_max'][m_S] = lofarcat['LR'][m_S]
    lofarcat['Ng_LR_good'][m_S] = 1*(lofarcat['LR_threshold'][m_S])

    # some flags for mult_gaus sources:
    # source has good LR match, and no gaus
    # multiple matches to different sources
    # source has no good LR match, but one gaus does
        
        
        
    lofarcat.write(lofarcat_file, overwrite=True, serialize_meta=True)
        
