import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, join, Column
import astropy.units as u
import astropy.coordinates as ac
import utils.plot_util as pp
import os

from lofar_source_sorter_dr2 import Mask


'''
M_Diagnosis_Code
 1 - accept ML source
 2 - accept ML gaus (there should be only one or best one...
 3 -
 4 - 
 5 - lgz
 6 -
>> 0: no match
>> 1: accept ML of the source
>> 2: accept ML of the gaussian with highest ML
>> 3: deblend and accept both gaussians
>> 4: deblend workflow
>> 5: LOFAR galaxy zoo
 
'''

if __name__=='__main__':
    
    size_large = 15.           # in arcsec
    separation1 = 45.          # in arcsec
    size_huge = 25.            # in arcsec
    #separation2 = 30.          # in arcsec
    fluxcut = 10               # in mJy
    fluxcut2 = 2.5               # in mJy


    step = int(sys.argv[1])
    if step not in  [1,2]:
        print('unknown step',step)
        sys.exit(1)

    mode = sys.argv[2]
    if mode.lower() not in ['isol', 'nonisol']:
        print('unknown mode',mode)
        sys.exit(1)
    #mode = 'isol'
    #mode = 'nonisol'

    path = '/data2/wwilliams/projects/lofar_surveys/LoTSS-DR2-Feb2020/'
    lofargcat_file = path+'LoTSS_DR2_rolling.gaus_0h.fits'
    lofarcat_file_srt = path+'LoTSS_DR2_rolling.srl_0h.sorted.fits'


    version = 'v1'
    
    if version == 'v1':
        lLR_thresh = 0.404            # LR threshold
        #lofarcat_file_srt = path+'LOFAR_HBA_T1_DR1_catalog_v0.99.srl.gmasked.sorted.v1.fits'
    else:
        print('Unknown version, quitting')
        sys.exit(1)


    #lofar_msource_flowchart_file = path + 'msources/LOFAR_flowchart.fixed.fits'
    if not os.path.exists(path + 'msources/'):
        os.mkdir(path + 'msources/')
    lofar_msource_flowchart_file = path + 'msources/step{st:d}_{m}_msources_flowchart_{v}.fits'.format(st=step, m=mode,v=version)

    # Gaus catalogue
    lofargcat = Table.read(lofargcat_file)
    # only relevant gaussians are in M or C sources
    lofargcat = lofargcat[lofargcat['S_Code'] != 'S']

    # Source catalogue
    lofarcat_full = Table.read(lofarcat_file_srt)

    # from Lara/Philip - classifications of compact isolated m sources
    #lofar_msource_flowchart = Table.read(lofar_msource_flowchart_file)

    if mode == 'isol':
        FFlag = 6
    elif mode == 'nonisol':
        FFlag = 8

    ## select only relevant sources - from step 1 of the flowchart (rerunning in step2 changes final outputs...)
    lofarcat = lofarcat_full[lofarcat_full['FC_flag1'] == FFlag]


    lofarcat.add_column(Column(99*np.ones(len(lofarcat),dtype=int), 'M_Diagnosis_Code'))

    source_nlr = lofarcat['LR'] < lLR_thresh
    source_lr = lofarcat['LR'] >= lLR_thresh
    source_lr2 = lofarcat['LR'] >= 10*lLR_thresh
    
    source_nglr = lofarcat['Ng_LR_good'] == 0
    source_glr = lofarcat['Ng_LR_good'] > 0
    source_glr1 = lofarcat['Ng_LR_good'] == 1
    source_glr2 = lofarcat['Ng_LR_good'] == 2
    source_glr3p = lofarcat['Ng_LR_good'] >= 3
    
    glarge_sep = lofarcat['G_max_sep'] >= 15.
    
    nncomplex = (lofarcat['NN_sep'] < 100.) & (lofarcat['NN_Maj'] > 10.)
    nnartefact = (lofarcat_full['Artefact_flag'][lofarcat['NN_idx']])
    
    m_S_small = (lofarcat['S_Code'] !='S') & (lofarcat['Maj'] <= 15.)
    m_S_small_isol = m_S_small & (lofarcat['NN_sep'] > 45.)
    
    #(lofarcat['Ng_LR_good'][m_S_small_isol] == 0 ) & (source_nlr[m_S_small_isol])
    
    #m_S_small_isol & source_nlr
    
    
    #m_S_small_isol & source_lr
    
    if step == 1:
        stvcol = 'green'
    elif step == 2:
        stvcol = 'black'
    
    masterlist = []
    
    M_small_m = Mask(lofarcat['FC_flag1'] == FFlag,
                    'small_m',
                    qlabel='small\nnon-isolated\nm\n',
                    masterlist=masterlist)
    
    
    #lofarcat['NN_sep'] > 45.
    M_small_m_nisol = M_small_m.submask(lofarcat['FC_flag1'] == FFlag,
                                        'nisol',
                                        qlabel='LRs >= 0.404??',
                                        edgelabel='Y',
                                        masterlist=masterlist)
    
    M_small_m_nisol_lr = M_small_m_nisol.submask(source_lr,
                                        'lr',
                                        edgelabel='Y (SLR)',
                                        qlabel='N(LRg >= 0.404) >= 1?',
                                        masterlist=masterlist)
    
    
    M_small_m_nisol_lr_glr = M_small_m_nisol_lr.submask(source_glr,
                                        'glr',
                                        edgelabel='Y (GLRSLR)',
                                        qlabel='N guassian lr?',
                                        masterlist=masterlist)
    
    
    M_small_m_nisol_lr_glr_glr1 = M_small_m_nisol_lr_glr.submask(source_glr1,
                                        'glr1',
                                        edgelabel='1',
                                        qlabel='same as source?',
                                        masterlist=masterlist)
    
    M_small_m_nisol_lr_glr_glr1_ges = M_small_m_nisol_lr_glr_glr1.submask(lofarcat['N_G_LR_matchsource'] == 1,
                                        'ges',
                                        edgelabel='Y',
                                        qlabel='accept ML source',
                                        color='blue',
                                        masterlist=masterlist)
    lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr1_ges.mask] = 1
    
    M_small_m_nisol_lr_glr_glr1_nges = M_small_m_nisol_lr_glr_glr1.submask(lofarcat['N_G_LR_matchsource'] != 1,
                                        'nges',
                                        edgelabel='N',
                                        qlabel='compare LRs and LRg',
                                        masterlist=masterlist)
    
    
    M_small_m_nisol_lr_glr_glr1_nges_1 = M_small_m_nisol_lr_glr_glr1_nges.submask(lofarcat['G_LR_case1'] == 1,
                                        'case1',
                                        edgelabel='1',
                                        qlabel='accept ML source',
                                        color='blue',
                                        masterlist=masterlist)
    lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr1_nges_1.mask] = 1
    
    
    M_small_m_nisol_lr_glr_glr1_nges_2 = M_small_m_nisol_lr_glr_glr1_nges.submask(lofarcat['G_LR_case1'] == 2,
                                        'case2',
                                        edgelabel='2',
                                        qlabel='accept ML gaus',
                                        color='blue',
                                        masterlist=masterlist)
    lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr1_nges_2.mask] = 2
    
    
    M_small_m_nisol_lr_glr_glr1_nges_3 = M_small_m_nisol_lr_glr_glr1_nges.submask(lofarcat['G_LR_case1'] == 3,
                                        'case3',
                                        edgelabel='3',
                                        qlabel='visual check',
                                        color=stvcol,
                                        masterlist=masterlist)
    
    # lofarcat['m_nisol_flag_vc2'] = {1,4,5} = lgz, match, deblend
    
    if step == 2:
        M_small_m_nisol_lr_glr_glr1_nges_3_lgz = M_small_m_nisol_lr_glr_glr1_nges_3.submask(lofarcat['m_nisol_flag_vc2'] == 1,
                                            'lgz',
                                            edgelabel='lgz',
                                            qlabel='lgz',
                                            color='green',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr1_nges_3_lgz.mask] = 5
        
        
        M_small_m_nisol_lr_glr_glr1_nges_3_blend = M_small_m_nisol_lr_glr_glr1_nges_3.submask(lofarcat['m_nisol_flag_vc2'] == 5,
                                            'blend',
                                            edgelabel='blend',
                                            qlabel='blend',
                                            color='orange',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr1_nges_3_blend.mask] = 4
        
        
        M_small_m_nisol_lr_glr_glr1_nges_3_nmatch = M_small_m_nisol_lr_glr_glr1_nges_3.submask(lofarcat['m_nisol_flag_vc2'] == 3,
                                            'match',
                                            edgelabel='no match',
                                            qlabel='accept no match',
                                            color='red',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr1_nges_3_nmatch.mask] = 0
        
        
        M_small_m_nisol_lr_glr_glr1_nges_3_match = M_small_m_nisol_lr_glr_glr1_nges_3.submask(lofarcat['m_nisol_flag_vc2'] == 4,
                                            'match',
                                            edgelabel='match',
                                            qlabel='accept match',
                                            color='blue',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr1_nges_3_match.mask] = 1
        
        
        M_small_m_nisol_lr_glr_glr1_nges_3_prob = M_small_m_nisol_lr_glr_glr1_nges_3.submask(lofarcat['m_nisol_flag_vc2'] == 0,
                                            'prob',
                                            edgelabel='prob',
                                            qlabel='prob',
                                            #color='',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr1_nges_3_prob.mask] = 6
    
    
    M_small_m_nisol_lr_glr_glr3p = M_small_m_nisol_lr_glr.submask(source_glr3p,
                                        'glr3p',
                                        qlabel='all G match S?',
                                        edgelabel='3+',
                                        masterlist=masterlist)

    
    M_small_m_nisol_lr_glr_glr3p_ges = M_small_m_nisol_lr_glr_glr3p.submask(lofarcat['N_G_LR_matchsource'] == lofarcat['Ng_LR_good'],
                                        'ges',
                                        edgelabel='Y',
                                        color='blue',
                                        qlabel='accept ML source',
                                        masterlist=masterlist)
    lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr3p_ges.mask] = 1
    
    M_small_m_nisol_lr_glr_glr3p_nges = M_small_m_nisol_lr_glr_glr3p.submask(lofarcat['N_G_LR_matchsource'] != lofarcat['Ng_LR_good'],
                                        'nges',
                                        edgelabel='N',
                                        qlabel='deblend',
                                        color='orange',
                                        masterlist=masterlist)
    lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr3p_nges.mask] = 4
    
    
    M_small_m_nisol_lr_glr_glr2 = M_small_m_nisol_lr_glr.submask(source_glr2,
                                        'glr2',
                                        'compare G/S matches',
                                        edgelabel='2',
                                        masterlist=masterlist)
    
    
    M_small_m_nisol_lr_glr_glr2_ges = M_small_m_nisol_lr_glr_glr2.submask(lofarcat['N_G_LR_matchsource'] == 2,
                                        'ges',
                                        edgelabel='both G match S',
                                        qlabel='accept ML source',
                                        color='blue',
                                        masterlist=masterlist)
    lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr2_ges.mask] = 1
    
    M_small_m_nisol_lr_glr_glr2_nges = M_small_m_nisol_lr_glr_glr2.submask(lofarcat['N_G_LR_matchsource'] == 0,
                                        'nges',
                                        edgelabel='neither G match S',
                                        qlabel='visual check',
                                        color=stvcol,
                                        masterlist=masterlist)
    
    # lofarcat['m_nisol_flag_vc2'] = {1,4,5} = lgz, match, deblend
    if step == 2:
        M_small_m_nisol_lr_glr_glr2_nges_lgz = M_small_m_nisol_lr_glr_glr2_nges.submask(lofarcat['m_nisol_flag_vc2'] == 1,
                                            'lgz',
                                            edgelabel='lgz',
                                            qlabel='lgz',
                                            color='green',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr2_nges_lgz.mask] = 5
        
        M_small_m_nisol_lr_glr_glr2_nges_blend = M_small_m_nisol_lr_glr_glr2_nges.submask(lofarcat['m_nisol_flag_vc2'] == 5,
                                            'blend',
                                            edgelabel='blend',
                                            qlabel='blend',
                                            color='orange',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr2_nges_blend.mask] = 4
        
        M_small_m_nisol_lr_glr_glr2_nges_match = M_small_m_nisol_lr_glr_glr2_nges.submask(lofarcat['m_nisol_flag_vc2'] == 4,
                                            'match',
                                            edgelabel='match',
                                            qlabel='accept match',
                                            color='blue',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr2_nges_match.mask] = 1
        
        
        M_small_m_nisol_lr_glr_glr2_nges_prob = M_small_m_nisol_lr_glr_glr2_nges.submask(lofarcat['m_nisol_flag_vc2'] == 0,
                                            'prob',
                                            edgelabel='prob',
                                            qlabel='prob',
                                            #color='blue',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr2_nges_prob.mask] = 6
    
    
    M_small_m_nisol_lr_glr_glr2_nges1 = M_small_m_nisol_lr_glr_glr2.submask(lofarcat['N_G_LR_matchsource'] == 1,
                                        'nges1',
                                        edgelabel='One G matches S',
                                        qlabel='compare LRsource and LRgaus\'s',
                                        #color='red',
                                        masterlist=masterlist)
    
    M_small_m_nisol_lr_glr_glr2_nges1_case1 = M_small_m_nisol_lr_glr_glr2_nges1.submask(lofarcat['G_LR_case2'] == 1,
                                        'case1',
                                        edgelabel='1',
                                        qlabel='accept ML source',
                                        color='blue',
                                        masterlist=masterlist)
    lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr2_nges1_case1.mask] = 1
    
    
    M_small_m_nisol_lr_glr_glr2_nges1_case2 = M_small_m_nisol_lr_glr_glr2_nges1.submask(lofarcat['G_LR_case2'] == 2,
                                        'case2',
                                        edgelabel='2',
                                        qlabel='deblend (direct)',
                                        color='orange',
                                        masterlist=masterlist)
    lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr2_nges1_case2.mask] = 3
    
    
    M_small_m_nisol_lr_glr_glr2_nges1_case3 = M_small_m_nisol_lr_glr_glr2_nges1.submask(lofarcat['G_LR_case2'] == 3,
                                        'case3',
                                        edgelabel='3',
                                        qlabel='deblend',
                                        color='orange',
                                        masterlist=masterlist)
    lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_glr_glr2_nges1_case3.mask] = 4
    
    
    M_small_m_nisol_lr_nglr = M_small_m_nisol_lr.submask(source_nglr,
                                        'nglr',
                                        edgelabel='N (nGLRSLR)',
                                        qlabel='LRs >= 10*0.404?',
                                        masterlist=masterlist)
    
    M_small_m_nisol_lr_nglr_lr2 = M_small_m_nisol_lr_nglr.submask(source_lr2,
                                        'slr',
                                        edgelabel='Y',
                                        qlabel='accept ML source',
                                        color='blue',
                                        masterlist=masterlist)
    lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_nglr_lr2.mask] = 1
    
    M_small_m_nisol_lr_nglr_nlr2 = M_small_m_nisol_lr_nglr.submask(~source_lr2,
                                        'nslr',
                                        edgelabel='N',
                                        qlabel='visual check',
                                        color=stvcol,
                                        masterlist=masterlist)
    
   
    # lofarcat['m_nisol_flag_vc2'] = {1,3,4,6} = lgz, no match, match, artefact
    
    if step == 2:
        M_small_m_nisol_lr_nglr_nlr2_lgz = M_small_m_nisol_lr_nglr_nlr2.submask(lofarcat['m_nisol_flag_vc2'] == 1,
                                            'lgz',
                                            edgelabel='lgz',
                                            qlabel='lgz',
                                            color='green',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_nglr_nlr2_lgz.mask] = 5
        
        
        M_small_m_nisol_lr_nglr_nlr2_no_match = M_small_m_nisol_lr_nglr_nlr2.submask(lofarcat['m_nisol_flag_vc2'] == 3 ,
                                            'no_match',
                                            edgelabel='no_match',
                                            qlabel='accept no match',
                                            color='red',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_nglr_nlr2_no_match.mask] = 0
        
        M_small_m_nisol_lr_nglr_nlr2_match = M_small_m_nisol_lr_nglr_nlr2.submask(lofarcat['m_nisol_flag_vc2'] == 4,
                                            'match',
                                            edgelabel='match',
                                            qlabel='match',
                                            color='blue',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_nglr_nlr2_match.mask] = 1
        
        M_small_m_nisol_lr_nglr_nlr2_artefact = M_small_m_nisol_lr_nglr_nlr2.submask(lofarcat['m_nisol_flag_vc2'] == 6,
                                            'artefact',
                                            edgelabel='artefact',
                                            qlabel='artefact',
                                            color='grey',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_nglr_nlr2_artefact.mask] = -1
        
        M_small_m_nisol_lr_nglr_nlr2_prob = M_small_m_nisol_lr_nglr_nlr2.submask(lofarcat['m_nisol_flag_vc2'] == 0,
                                            'prob',
                                            edgelabel='prob',
                                            qlabel='prob',
                                            #color='grey',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_lr_nglr_nlr2_prob.mask] = 6
    
    
    
    
    
    M_small_m_nisol_nlr = M_small_m_nisol.submask(source_nlr,
                                        'nlr',
                                        edgelabel='N (nSLR)',
                                        qlabel='N(LRg >= 0.404) >= 1?',
                                        masterlist=masterlist)
    

    M_small_m_nisol_nlr_glr = M_small_m_nisol_nlr.submask(source_glr,
                                        'glr',
                                        edgelabel='Y (GLRnSLR)',
                                        qlabel='N(LRg >= 0.404)?',
                                        masterlist=masterlist)
    
    
    M_small_m_nisol_nlr_glr_glr1 = M_small_m_nisol_nlr_glr.submask(source_glr1,
                                        'glr1',
                                        edgelabel='1',
                                        qlabel='LRgaus > 10*0.404 & \nMaj_gaus< 10"?',
                                        masterlist=masterlist)
    
    M_small_m_nisol_nlr_glr_glr1_case1 = M_small_m_nisol_nlr_glr_glr1.submask(lofarcat['G_LR_case3']==1,
                                        'case1',
                                        edgelabel='Y',
                                        qlabel='accept ML gaus',
                                        color='blue',
                                        masterlist=masterlist)
    lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_glr_glr1_case1.mask] = 2
    
    M_small_m_nisol_nlr_glr_glr1_case2 = M_small_m_nisol_nlr_glr_glr1.submask(lofarcat['G_LR_case3']==2,
                                        'case2',
                                        edgelabel='N',
                                        qlabel='visual check',
                                        color=stvcol,
                                        masterlist=masterlist)
    
    
    # lofarcat['m_nisol_flag_vc2'] = {1,3,4,5,6} = lgz, no match, match, deblend, artefact
    if step == 2:
        M_small_m_nisol_nlr_glr_glr1_case2_lgz = M_small_m_nisol_nlr_glr_glr1_case2.submask(lofarcat['m_nisol_flag_vc2']==1,
                                            'lgz',
                                            edgelabel='lgz',
                                            qlabel='lgz',
                                            color='green',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_glr_glr1_case2_lgz.mask] = 5
        
        
        M_small_m_nisol_nlr_glr_glr1_case2_no_match = M_small_m_nisol_nlr_glr_glr1_case2.submask(lofarcat['m_nisol_flag_vc2']==3,
                                            'no_match',
                                            edgelabel='no_match',
                                            qlabel='accept no match',
                                            color='red',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_glr_glr1_case2_no_match.mask] = 0
        
        M_small_m_nisol_nlr_glr_glr1_case2_match = M_small_m_nisol_nlr_glr_glr1_case2.submask(lofarcat['m_nisol_flag_vc2']==4,
                                            'match',
                                            edgelabel='match',
                                            qlabel='match',
                                            color='blue',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_glr_glr1_case2_match.mask] = 2
        
        M_small_m_nisol_nlr_glr_glr1_case2_deblend = M_small_m_nisol_nlr_glr_glr1_case2.submask(lofarcat['m_nisol_flag_vc2']==5,
                                            'deblend',
                                            edgelabel='deblend',
                                            qlabel='deblend',
                                            color='orange',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_glr_glr1_case2_deblend.mask] = 4
        
        M_small_m_nisol_nlr_glr_glr1_case2_artefact = M_small_m_nisol_nlr_glr_glr1_case2.submask(lofarcat['m_nisol_flag_vc2']==6,
                                            'artefact',
                                            edgelabel='artefact',
                                            qlabel='artefact',
                                            color='grey',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_glr_glr1_case2_artefact.mask] = -1
        
        
        
        M_small_m_nisol_nlr_glr_glr1_case2_prob = M_small_m_nisol_nlr_glr_glr1_case2.submask(lofarcat['m_nisol_flag_vc2']==0,
                                            'prob',
                                            edgelabel='prob',
                                            qlabel='prob',
                                            #color='grey',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_glr_glr1_case2_prob.mask] = 6
    
    
    
    M_small_m_nisol_nlr_glr_glr2 = M_small_m_nisol_nlr_glr.submask(source_glr2,
                                        'glr2',
                                        qlabel='match same gal?',
                                        edgelabel='2',
                                        masterlist=masterlist)
    
    
    M_small_m_nisol_nlr_glr_glr2_same = M_small_m_nisol_nlr_glr_glr2.submask(lofarcat['Ng_LR_good_unique'] == 1,
                                        'same',
                                        qlabel='visual check',
                                        color=stvcol,
                                        edgelabel='Y',
                                        masterlist=masterlist)
    
    
    # lofarcat['m_nisol_flag_vc2'] = {1,3,4} = lgz, no match, match
    
    
    if step == 2:
        M_small_m_nisol_nlr_glr_glr2_same_lgz = M_small_m_nisol_nlr_glr_glr2_same.submask(lofarcat['m_nisol_flag_vc2'] == 1,
                                            'lgz',
                                            qlabel='lgz',
                                            color='green',
                                            edgelabel='lgz',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_glr_glr2_same_lgz.mask] = 5
        
        M_small_m_nisol_nlr_glr_glr2_same_no_match = M_small_m_nisol_nlr_glr_glr2_same.submask(lofarcat['m_nisol_flag_vc2'] == 3,
                                            'no_match',
                                            qlabel='accept no match',
                                            color='red',
                                            edgelabel='no_match',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_glr_glr2_same_no_match.mask] = 0
        
        M_small_m_nisol_nlr_glr_glr2_same_match = M_small_m_nisol_nlr_glr_glr2_same.submask(lofarcat['m_nisol_flag_vc2'] == 4,
                                            'match',
                                            qlabel='match',
                                            color='blue',
                                            edgelabel='match',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_glr_glr2_same_match.mask] = 2
        
        
        M_small_m_nisol_nlr_glr_glr2_same_prob = M_small_m_nisol_nlr_glr_glr2_same.submask(lofarcat['m_nisol_flag_vc2'] == 0,
                                            'prob',
                                            qlabel='prob',
                                            #color='blue',
                                            edgelabel='prob',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_glr_glr2_same_prob.mask] = 6
    
    
    M_small_m_nisol_nlr_glr_glr2_diff = M_small_m_nisol_nlr_glr_glr2.submask(lofarcat['Ng_LR_good_unique'] != 1,
                                        'diff',
                                        qlabel='deblend',
                                        color='orange',
                                        edgelabel='N',
                                        masterlist=masterlist)
    lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_glr_glr2_diff.mask] = 4
    
    M_small_m_nisol_nlr_glr_glr3p = M_small_m_nisol_nlr_glr.submask(source_glr3p,
                                        'glr3p',
                                        edgelabel='3+',
                                        qlabel='lgz',
                                        color='green',
                                        masterlist=masterlist)
    lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_glr_glr3p.mask] = 5


    M_small_m_nisol_nlr_nglr = M_small_m_nisol_nlr.submask(source_nglr,
                                        'nglr',
                                        edgelabel='N   (nGLRnSLR)',
                                        qlabel='max(sep_g) < 15"?',
                                        masterlist=masterlist)
    
    
    
    M_small_m_nisol_nlr_nglr_nglargesep = M_small_m_nisol_nlr_nglr.submask(~glarge_sep,
                                        'nglargesep',
                                        edgelabel='Y',
                                        qlabel='NNsep<100" & NNsize>10"?',
                                        masterlist=masterlist)
    
    
    M_small_m_nisol_nlr_nglr_nglargesep_nnncomplex = M_small_m_nisol_nlr_nglr_nglargesep.submask(~nncomplex,
                                        'nnncomplex',
                                        edgelabel='N',
                                        qlabel='accept no match',
                                        color='red',
                                        masterlist=masterlist)
    lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_nglr_nglargesep_nnncomplex.mask] = 0
    
    
    M_small_m_nisol_nlr_nglr_nglargesep_nncomplex = M_small_m_nisol_nlr_nglr_nglargesep.submask(nncomplex,
                                        'nncomplex',
                                        edgelabel='Y',
                                        qlabel='NN is artefact?',
                                        masterlist=masterlist)
    
    
    M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_art = M_small_m_nisol_nlr_nglr_nglargesep_nncomplex.submask(nnartefact,
                                        'nnartefact',
                                        edgelabel='Y',
                                        qlabel='accept no match',
                                        color='red',
                                        masterlist=masterlist)
    lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_art.mask] = 0
    
    
    M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart = M_small_m_nisol_nlr_nglr_nglargesep_nncomplex.submask(~nnartefact,
                                        'nnnartefact',
                                        edgelabel='N',
                                        qlabel='visual check',
                                        color=stvcol,
                                        #color='green',
                                        masterlist=masterlist)
    #1 Complex (lgz)
    #2 complex (lgz - should be done already)
    #3 No match
    #4 artefact
    #5 no match, but NN is artefact
    #6 other (redo)
    #m_nisol_flag_vc1
    if step == 2:
        M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart_lgz = M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart.submask(lofarcat['m_nisol_flag_vc1']==1,
                                            'lgz',
                                            edgelabel='Y',
                                            qlabel='lgz',
                                            color='green',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart_lgz.mask] = 5 # this was 51
        
        M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart_lgz1 = M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart.submask(lofarcat['m_nisol_flag_vc1']==2,
                                            'lgz1',
                                            edgelabel='Y',
                                            qlabel='lgz (already)',
                                            color='green',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart_lgz1.mask] = 0  # most of these get associated...
        
        M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart_nlgz = M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart.submask(lofarcat['m_nisol_flag_vc1']>2,
                                            'nlgz',
                                            edgelabel='N',
                                            qlabel='no match',
                                            color='red',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart_nlgz.mask] = 0
        
        
        M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart_prob = M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart.submask(lofarcat['m_nisol_flag_vc1']==0,
                                            'prob',
                                            edgelabel='prob',
                                            qlabel='prob',
                                            #color='red',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart_prob.mask] = 6
    
    M_small_m_nisol_nlr_nglr_glargesep = M_small_m_nisol_nlr_nglr.submask(glarge_sep,
                                        'glargesep',
                                        edgelabel='N',
                                        qlabel='visual check',
                                        color=stvcol,
                                        masterlist=masterlist)
    
    # lofarcat['m_nisol_flag_vc2'] = {1,3,6} = lgz, no match, artefact
    
    # we've not done visual classification
    if step == 2:
        M_small_m_nisol_nlr_nglr_glargesep_lgz = M_small_m_nisol_nlr_nglr_glargesep.submask(lofarcat['m_nisol_flag_vc2']==1,
                                            'lgz',
                                            edgelabel='lgz',
                                            qlabel='lgz',
                                            color='green',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_nglr_glargesep_lgz.mask] = 5
        
        
        M_small_m_nisol_nlr_nglr_glargesep_no_match = M_small_m_nisol_nlr_nglr_glargesep.submask(lofarcat['m_nisol_flag_vc2']==3,
                                            'no_match',
                                            edgelabel='no_match',
                                            qlabel='accept no match',
                                            color='red',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_nglr_glargesep_no_match.mask] = 0
        
        
        M_small_m_nisol_nlr_nglr_glargesep_artefact = M_small_m_nisol_nlr_nglr_glargesep.submask(lofarcat['m_nisol_flag_vc2']==6,
                                            'artefact',
                                            edgelabel='artefact',
                                            qlabel='artefact',
                                            color='grey',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_nglr_glargesep_artefact.mask] = -1
        
        M_small_m_nisol_nlr_nglr_glargesep_prob = M_small_m_nisol_nlr_nglr_glargesep.submask(lofarcat['m_nisol_flag_vc2']==0,
                                            'prob',
                                            edgelabel='prob',
                                            qlabel='prob',
                                            #color='grey',
                                            masterlist=masterlist)
        lofarcat['M_Diagnosis_Code'][M_small_m_nisol_nlr_nglr_glargesep_prob.mask] = 6
        
        
    #if 'FC_flag1' not in lofarcat.colnames:
        #lofarcat.add_column(Column(-1*np.ones(len(lofarcat),dtype=int), 'FC_flag1'))
    #else:
        #lofarcat['FC_flag1']=Column(-1*np.ones(len(lofarcat),dtype=int))



    mcflg = 'MC_flag{st:d}'.format(st=step)
    if mcflg not in lofarcat.colnames:
        lofarcat.add_column(Column(-1*np.ones(len(lofarcat),dtype=int), mcflg))
    else:
        lofarcat[mcflg] = -1*np.ones(len(lofarcat),dtype=int)
    i = 0
    for t in masterlist:
        if not t.has_children:
            lofarcat[mcflg][t.mask] = i
            i += 1
            
    
    # make flowchart from list of masks
    plot_flowchart = True
    plot_verbose = False
    try:
        import pygraphviz as pgv
    except ImportError:
        print('no pygraphviz; cannot make visual flowchart')
        plot_flowchart = False
    if plot_flowchart:

        PW = 75.

        A=pgv.AGraph(directed=True, strict=True)
        A.edge_attr['arrowhead']='none'
        A.node_attr['start']='south'
        A.node_attr['end']='north'
        A.node_attr['style']='filled'
        A.node_attr['fillcolor']='white'
        A.node_attr['fixed-size']='true'
        A.node_attr['width']='2.5'
        A.node_attr['height']='2'
        A.edge_attr['color']='gray'
        A.edge_attr['tailclip']='false'
        A.edge_attr['headclip']='false'
        A.graph_attr['outputorder'] = 'edgesfirst'
        #A.graph_attr['splines'] = 'ortho'  # orthogonal
        A.graph_attr['rankdir'] = 'TB'

        #A.add_node('start', label='ALL\n{n:n}'.format(n=M_small_m.N), shape='parallelogram') 
        #A.add_node('m_all', label='Large?'.format(n=M_small_m.N), shape='diamond') 

        #A.add_edge('start', 'all', label='', penwidth=M_small_m.f*PW)
        i = 0
        for t in masterlist:
            
            
            if t.has_children:
                shape='diamond'         # intermediate point is a question
                
                if t.p < 1.:
                    label='{lab:s}\n{n:n}\n{p:.2f}%'.format(lab=t.qlabel,n=t.n,p=t.p)
                if t.p < 10.:
                    label='{lab:s}\n{n:n}\n{p:.1f}%'.format(lab=t.qlabel,n=t.n,p=t.p)
                else:
                    label='{lab:s}\n{n:n}\n{p:.0f}%'.format(lab=t.qlabel,n=t.n,p=t.p)
            else:
                shape='parallelogram'   # end point is a final mask
                if t.p < 1.:
                    label='- {i:n} -\n{lab:s}\n{n:n}\n{p:.2f}%'.format(i=i,lab=t.qlabel,n=t.n,p=t.p)
                elif t.p < 10.:
                    label='- {i:n} -\n{lab:s}\n{n:n}\n{p:.1f}%'.format(i=i,lab=t.qlabel,n=t.n,p=t.p)
                else:
                    label='- {i:n} -\n{lab:s}\n{n:n}\n{p:.0f}%'.format(i=i,lab=t.qlabel,n=t.n,p=t.p)
                i += 1
            if t.color:
                c = t.color
            else:
                c = 'black'
            # add node
            A.add_node(t.name, label=label, shape=shape, color=c)
            
            # add edge to parent
            if t.has_parent:
                A.add_edge(t.parent.name, t.name, label=t.edgelabel, penwidth=t.f*PW)

        if plot_verbose:
            print((A.string())) # print dot file to standard output

        # make the flowchart
        #Optional prog=['neato'|'dot'|'twopi'|'circo'|'fdp'|'nop']
        #neato, dot, twopi, circo, fdp, nop, wc, acyclic, gvpr, gvcolor, ccomps, sccmap, tred, sfdp.
        A.layout('dot') # layout with dot
        A.draw('flow_step{st:d}_s{s:.0f}_nn{nn:.0f}_msources_{m}_{v}.png'.format(st=step,s=size_large,nn=separation1,v=version,m=mode)) # write to file
        A.draw('flow_step{st:d}_s{s:.0f}_nn{nn:.0f}_msources_{m}_{v}.pdf'.format(st=step,s=size_large,nn=separation1,v=version,m=mode)) # write to file
        A.write('flow_step{st:d}_s{s:.0f}_nn{nn:.0f}_msources_{m}_{v}.dot'.format(st=step,s=size_large,nn=separation1,v=version,m=mode)) # write to file

        
    # make a test sample for each final mask
    makesample = 1
    if makesample:
        for t in masterlist:
            if not t.has_children :
                print(t.name)
                t.make_sample(lofarcat,Nsample=None)
        
        
    lofarcat.keep_columns(['Source_Name','M_Diagnosis_Code',mcflg])
    lofarcat[M_small_m_nisol.mask].write(lofar_msource_flowchart_file, overwrite=True)
