import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from astropy.table import Table, join, Column
import astropy.units as u
import astropy.coordinates as ac
#import utils.plot_util as pp
import utils.plotting as pp

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
>> 6: visual check (early steps)
 
'''


def update_flags(lofarcat_file_srt, msource_cat, stepi):
    ''' code originally in get_msource_flags '''
    

    lofarcat = Table.read(lofarcat_file_srt)



    #################################################################################
    # the following come from outputs from Lara/Philip for compact isolated m sources

    #################################################################################
    ### msource1_flag
    #0: no match
    #1: accept ML of the source
    #2: accept ML of the gaussian with highest ML
    #3: deblend and accept both gaussians
    #4: deblend workflow
    #5: LOFAR galaxy zoo 
    #5: visual check 


    print ('saving flags to lofar cat')
    #if 'msource_flag{:d}'.format(stepi) in lofarcat.colnames:
        #lofarcat.remove_column('msource_flag{:d}'.format(flagi))
    #if 'MC_flag{:d}'.format(stepi) in lofarcat.colnames:
        #lofarcat.remove_column('MC_flag{:d}'.format(stepi))
    lofarcat.sort('Source_Name')
    tt=join(lofarcat, msource_cat, join_type='left', keys=['Source_Name'])
    if 'msource_flag{st:d}_2'.format(st=stepi) in tt.colnames:
        tt.rename_column('msource_flag{st:d}_2'.format(st=stepi), 'msource_flag{st:d}'.format(st=stepi))
        tt.rename_column('MC_flag{st:d}_2'.format(st=stepi), 'MC_flag{st:d}'.format(st=stepi))
        
    tt['msource_flag{st:d}'.format(st=stepi)].fill_value = -1
    tt['MC_flag{st:d}'.format(st=stepi)].fill_value = -1
    tt = tt.filled()
    tt.sort('Source_Name')

    if 'msource_flag{:d}'.format(stepi) in lofarcat.colnames:
        indm = (tt['MC_flag{:d}'.format(stepi)] > -1)
        lofarcat['MC_flag{:d}'.format(stepi)][indm] =  tt['MC_flag{:d}'.format(stepi)][indm]
        lofarcat['msource_flag{:d}'.format(stepi)][indm] =  tt['msource_flag{:d}'.format(stepi)][indm]

    else:
        lofarcat.add_column(tt['msource_flag{:d}'.format(stepi)])
        lofarcat.add_column(tt['MC_flag{:d}'.format(stepi)])


    #################################################################################

    ## write output file
    lofarcat.write(lofarcat_file_srt, overwrite=True, serialize_meta=True)
    
    
    return


if __name__=='__main__':
    
    size_large = 15.           # in arcsec
    separation1 = 45.          # in arcsec
    size_huge = 25.            # in arcsec
    
    #fluxcut2 = 4.0               # in mJy
    fluxcut2 = 4.  #float(sys.argv[4])
    #fluxcut2 = 3.0               # in mJy
    #fluxcut2 = 5.0               # in mJy
    #separation2 = 30.          # in arcsec

    if len(sys.argv) == 1:
        print("Usage is : python handle_m_sources_dr2.py field_code step_no mode")
        print('E.g.: python handle_m_sources_dr2.py 0 1 nonisol')
        sys.exit(1)

    h = str(sys.argv[1])
    if 'h' not in h:
        h+='h'
    if h not in  ['0h','13h']:
        print('unknown field code (should be 0h or 13h)',h)
        sys.exit(1)

    step = int(sys.argv[2])
    if step not in  [1,2]:
        print('unknown step',step)
        sys.exit(1)


    modes = ['isol', 'nonisol', 'isol-faint', 'nonisol-faint']

    mode = sys.argv[3]
    if mode.lower() not in ['all'] + modes:
        print('unknown mode',mode)
        sys.exit(1)
    #mode = 'isol'
    #mode = 'nonisol'
    if mode.lower() != 'all':
        modes = [mode]

    version ='v100'
    
    path = '/data2/wwilliams/projects/lofar_surveys/LoTSS-DR2-Feb2020/'
    #lofargcat_file = path+'lr/LoTSS_DR2_{version}.gaus_{h}.lr-full.fits'.format(h=h,version=version)
    lofarcat_file_srt = path+'LoTSS_DR2_{version}.srl_{h}.lr-full.sorted_step1_flux{ff:.0f}.hdf5'.format(h=h,version=version,ff=fluxcut2)

    if not os.path.exists(lofarcat_file_srt):
        print('{f} does not exist - first run lofar_source_sorter_dr2.py to produce this'.format(f=lofarcat_file_srt))

    

    for mode in modes:

        #lofar_msource_flowchart_file = path + 'msources/LOFAR_flowchart.fixed.fits'
        if not os.path.exists(path + 'msources/'):
            os.mkdir(path + 'msources/')
        lofar_msource_flowchart_file = path + 'msources/{h}_step{st:d}_{m}_msources_flowchart_{v}.fits'.format(h=h,st=step, m=mode,v=version)

        ## Gaus catalogue
        #lofargcat = Table.read(lofargcat_file, memmap=True)
        ## only relevant gaussians are in M or C sources
        #lofargcat = lofargcat[lofargcat['S_Code'] != 'S']

        # Source catalogue
        lofarcat_full = Table.read(lofarcat_file_srt)

        # from Lara/Philip - classifications of compact isolated m sources
        #lofar_msource_flowchart = Table.read(lofar_msource_flowchart_file)

        # NOTE : update based on latest main flowchart - note the difference for the two fields means different FC_flags
        if h=='0h':
            if mode == 'isol':
                FFlag = 10
            elif mode == 'isol-faint':
                FFlag = 11
            elif mode == 'nonisol':
                FFlag = 14
            elif mode == 'nonisol-faint':
                FFlag = 15
        elif h == '13h':
            if mode == 'isol':
                FFlag = 9
            elif mode == 'isol-faint':
                FFlag = 10
            elif mode == 'nonisol':
                FFlag = 14
            elif mode == 'nonisol-faint':
                FFlag = 15

        ## select only relevant sources - from step 1 of the flowchart (rerunning in step2 changes final outputs...)
        lofarcat = lofarcat_full[lofarcat_full['FC_flag1'] == FFlag]

        
        msourceflg = 'msource_flag{st:d}'.format(st=step)  # used to be M_Diagnosis_Code
        if msourceflg not in lofarcat.colnames:
            lofarcat.add_column(Column(99*np.ones(len(lofarcat),dtype=int), msourceflg))

        source_lr = lofarcat['LR_threshold']
        
        source_nlr = ~source_lr
        
        source_lr2 = lofarcat['LR_threshold10']

        
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
        
        # remove reference to full cat to save memory
        del lofarcat_full
        
        if step in [1,2]:
            stvcol = 'orange'
        elif step == 3:
            stvcol = 'black'
        
        masterlist = []
        
        M_small_m = Mask(lofarcat['FC_flag1'] == FFlag,
                        'small_m',
                        qlabel='small\n{m}\nm\n'.format(m=mode),
                        masterlist=masterlist,
                        maskflag=lofarcat['ML_flag'], maskflagstr='LR:LGZ')
        
        
        #lofarcat['NN_sep'] > 45.
        M_small_m_nisol = M_small_m.submask(lofarcat['FC_flag1'] == FFlag,
                                            'nisol',
                                            qlabel='LRs >= thresh??',
                                            edgelabel='Y',
                                            masterlist=masterlist)
        
        M_small_m_nisol_lr = M_small_m_nisol.submask(source_lr,
                                            'lr',
                                            edgelabel='Y (SLR)',
                                            qlabel='N(LRg >= thresh) >= 1?',
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
        lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr1_ges.mask] = 1
        
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
        lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr1_nges_1.mask] = 1
        
        
        M_small_m_nisol_lr_glr_glr1_nges_2 = M_small_m_nisol_lr_glr_glr1_nges.submask(lofarcat['G_LR_case1'] == 2,
                                            'case2',
                                            edgelabel='2',
                                            qlabel='accept ML gaus',
                                            color='blue',
                                            masterlist=masterlist)
        lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr1_nges_2.mask] = 2
        
        
        M_small_m_nisol_lr_glr_glr1_nges_3 = M_small_m_nisol_lr_glr_glr1_nges.submask(lofarcat['G_LR_case1'] == 3,
                                            'case3',
                                            edgelabel='3',
                                            qlabel='visual check',
                                            color=stvcol,
                                            masterlist=masterlist)
        lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr1_nges_3.mask] = 6
        
        # lofarcat['m_nisol_flag_vc2'] = {1,4,5} = lgz, match, deblend
        
        if step == 3:
            M_small_m_nisol_lr_glr_glr1_nges_3_lgz = M_small_m_nisol_lr_glr_glr1_nges_3.submask(lofarcat['m_nisol_flag_vc2'] == 1,
                                                'lgz',
                                                edgelabel='lgz',
                                                qlabel='lgz',
                                                color='orange',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr1_nges_3_lgz.mask] = 5
            
            
            M_small_m_nisol_lr_glr_glr1_nges_3_blend = M_small_m_nisol_lr_glr_glr1_nges_3.submask(lofarcat['m_nisol_flag_vc2'] == 5,
                                                'blend',
                                                edgelabel='blend',
                                                qlabel='blend',
                                                color='yellow',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr1_nges_3_blend.mask] = 4
            
            
            M_small_m_nisol_lr_glr_glr1_nges_3_nmatch = M_small_m_nisol_lr_glr_glr1_nges_3.submask(lofarcat['m_nisol_flag_vc2'] == 3,
                                                'match',
                                                edgelabel='no match',
                                                qlabel='accept no match',
                                                color='red',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr1_nges_3_nmatch.mask] = 0
            
            
            M_small_m_nisol_lr_glr_glr1_nges_3_match = M_small_m_nisol_lr_glr_glr1_nges_3.submask(lofarcat['m_nisol_flag_vc2'] == 4,
                                                'match',
                                                edgelabel='match',
                                                qlabel='accept match',
                                                color='blue',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr1_nges_3_match.mask] = 1
            
            
            M_small_m_nisol_lr_glr_glr1_nges_3_prob = M_small_m_nisol_lr_glr_glr1_nges_3.submask(lofarcat['m_nisol_flag_vc2'] == 0,
                                                'prob',
                                                edgelabel='prob',
                                                qlabel='prob',
                                                #color='',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr1_nges_3_prob.mask] = 6
        
        
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
        lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr3p_ges.mask] = 1
        
        M_small_m_nisol_lr_glr_glr3p_nges = M_small_m_nisol_lr_glr_glr3p.submask(lofarcat['N_G_LR_matchsource'] != lofarcat['Ng_LR_good'],
                                            'nges',
                                            edgelabel='N',
                                            qlabel='deblend',
                                            color='yellow',
                                            masterlist=masterlist)
        lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr3p_nges.mask] = 4
        
        
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
        lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr2_ges.mask] = 1
        
        M_small_m_nisol_lr_glr_glr2_nges = M_small_m_nisol_lr_glr_glr2.submask(lofarcat['N_G_LR_matchsource'] == 0,
                                            'nges',
                                            edgelabel='neither G match S',
                                            qlabel='visual check',
                                            color=stvcol,
                                            masterlist=masterlist)
        lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr2_nges.mask] = 6
        
        # lofarcat['m_nisol_flag_vc2'] = {1,4,5} = lgz, match, deblend
        if step == 3:
            M_small_m_nisol_lr_glr_glr2_nges_lgz = M_small_m_nisol_lr_glr_glr2_nges.submask(lofarcat['m_nisol_flag_vc2'] == 1,
                                                'lgz',
                                                edgelabel='lgz',
                                                qlabel='lgz',
                                                color='orange',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr2_nges_lgz.mask] = 5
            
            M_small_m_nisol_lr_glr_glr2_nges_blend = M_small_m_nisol_lr_glr_glr2_nges.submask(lofarcat['m_nisol_flag_vc2'] == 5,
                                                'blend',
                                                edgelabel='blend',
                                                qlabel='blend',
                                                color='yellow',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr2_nges_blend.mask] = 4
            
            M_small_m_nisol_lr_glr_glr2_nges_match = M_small_m_nisol_lr_glr_glr2_nges.submask(lofarcat['m_nisol_flag_vc2'] == 4,
                                                'match',
                                                edgelabel='match',
                                                qlabel='accept match',
                                                color='blue',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr2_nges_match.mask] = 1
            
            
            M_small_m_nisol_lr_glr_glr2_nges_prob = M_small_m_nisol_lr_glr_glr2_nges.submask(lofarcat['m_nisol_flag_vc2'] == 0,
                                                'prob',
                                                edgelabel='prob',
                                                qlabel='prob',
                                                #color='blue',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr2_nges_prob.mask] = 6
        
        
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
        lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr2_nges1_case1.mask] = 1
        
        
        M_small_m_nisol_lr_glr_glr2_nges1_case2 = M_small_m_nisol_lr_glr_glr2_nges1.submask(lofarcat['G_LR_case2'] == 2,
                                            'case2',
                                            edgelabel='2',
                                            qlabel='deblend (direct)',
                                            color='yellow',
                                            masterlist=masterlist)
        lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr2_nges1_case2.mask] = 3
        
        
        M_small_m_nisol_lr_glr_glr2_nges1_case3 = M_small_m_nisol_lr_glr_glr2_nges1.submask(lofarcat['G_LR_case2'] == 3,
                                            'case3',
                                            edgelabel='3',
                                            qlabel='deblend',
                                            color='yellow',
                                            masterlist=masterlist)
        lofarcat[msourceflg][M_small_m_nisol_lr_glr_glr2_nges1_case3.mask] = 4
        
        
        M_small_m_nisol_lr_nglr = M_small_m_nisol_lr.submask(source_nglr,
                                            'nglr',
                                            edgelabel='N (nGLRSLR)',
                                            qlabel='LRs >= 10*thresh?',
                                            masterlist=masterlist)
        
        M_small_m_nisol_lr_nglr_lr2 = M_small_m_nisol_lr_nglr.submask(source_lr2,
                                            'slr',
                                            edgelabel='Y',
                                            qlabel='accept ML source',
                                            color='blue',
                                            masterlist=masterlist)
        lofarcat[msourceflg][M_small_m_nisol_lr_nglr_lr2.mask] = 1
        
        M_small_m_nisol_lr_nglr_nlr2 = M_small_m_nisol_lr_nglr.submask(~source_lr2,
                                            'nslr',
                                            edgelabel='N',
                                            qlabel='visual check',
                                            color=stvcol,
                                            masterlist=masterlist)
        lofarcat[msourceflg][M_small_m_nisol_lr_nglr_nlr2.mask] = 6
        
    
        # lofarcat['m_nisol_flag_vc2'] = {1,3,4,6} = lgz, no match, match, artefact
        
        if step == 3:
            M_small_m_nisol_lr_nglr_nlr2_lgz = M_small_m_nisol_lr_nglr_nlr2.submask(lofarcat['m_nisol_flag_vc2'] == 1,
                                                'lgz',
                                                edgelabel='lgz',
                                                qlabel='lgz',
                                                color='orange',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_lr_nglr_nlr2_lgz.mask] = 5
            
            
            M_small_m_nisol_lr_nglr_nlr2_no_match = M_small_m_nisol_lr_nglr_nlr2.submask(lofarcat['m_nisol_flag_vc2'] == 3 ,
                                                'no_match',
                                                edgelabel='no_match',
                                                qlabel='accept no match',
                                                color='red',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_lr_nglr_nlr2_no_match.mask] = 0
            
            M_small_m_nisol_lr_nglr_nlr2_match = M_small_m_nisol_lr_nglr_nlr2.submask(lofarcat['m_nisol_flag_vc2'] == 4,
                                                'match',
                                                edgelabel='match',
                                                qlabel='match',
                                                color='blue',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_lr_nglr_nlr2_match.mask] = 1
            
            M_small_m_nisol_lr_nglr_nlr2_artefact = M_small_m_nisol_lr_nglr_nlr2.submask(lofarcat['m_nisol_flag_vc2'] == 6,
                                                'artefact',
                                                edgelabel='artefact',
                                                qlabel='artefact',
                                                color='grey',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_lr_nglr_nlr2_artefact.mask] = -1
            
            M_small_m_nisol_lr_nglr_nlr2_prob = M_small_m_nisol_lr_nglr_nlr2.submask(lofarcat['m_nisol_flag_vc2'] == 0,
                                                'prob',
                                                edgelabel='prob',
                                                qlabel='prob',
                                                #color='grey',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_lr_nglr_nlr2_prob.mask] = 6
        
        
        
        
        
        M_small_m_nisol_nlr = M_small_m_nisol.submask(source_nlr,
                                            'nlr',
                                            edgelabel='N (nSLR)',
                                            qlabel='N(LRg >= thresh) >= 1?',
                                            masterlist=masterlist)
        

        M_small_m_nisol_nlr_glr = M_small_m_nisol_nlr.submask(source_glr,
                                            'glr',
                                            edgelabel='Y (GLRnSLR)',
                                            qlabel='N(LRg >= thresh)?',
                                            masterlist=masterlist)
        
        
        M_small_m_nisol_nlr_glr_glr1 = M_small_m_nisol_nlr_glr.submask(source_glr1,
                                            'glr1',
                                            edgelabel='1',
                                            qlabel='LRgaus > 10*thresh & \nMaj_gaus< 10"?',
                                            masterlist=masterlist)
        
        M_small_m_nisol_nlr_glr_glr1_case1 = M_small_m_nisol_nlr_glr_glr1.submask(lofarcat['G_LR_case3']==1,
                                            'case1',
                                            edgelabel='Y',
                                            qlabel='accept ML gaus',
                                            color='blue',
                                            masterlist=masterlist)
        lofarcat[msourceflg][M_small_m_nisol_nlr_glr_glr1_case1.mask] = 2
        
        M_small_m_nisol_nlr_glr_glr1_case2 = M_small_m_nisol_nlr_glr_glr1.submask(lofarcat['G_LR_case3']==2,
                                            'case2',
                                            edgelabel='N',
                                            qlabel='visual check',
                                            color=stvcol,
                                            masterlist=masterlist)
        lofarcat[msourceflg][M_small_m_nisol_nlr_glr_glr1_case2.mask] = 6
        
        
        # lofarcat['m_nisol_flag_vc2'] = {1,3,4,5,6} = lgz, no match, match, deblend, artefact
        if step == 3:
            M_small_m_nisol_nlr_glr_glr1_case2_lgz = M_small_m_nisol_nlr_glr_glr1_case2.submask(lofarcat['m_nisol_flag_vc2']==1,
                                                'lgz',
                                                edgelabel='lgz',
                                                qlabel='lgz',
                                                color='orange',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_glr_glr1_case2_lgz.mask] = 5
            
            
            M_small_m_nisol_nlr_glr_glr1_case2_no_match = M_small_m_nisol_nlr_glr_glr1_case2.submask(lofarcat['m_nisol_flag_vc2']==3,
                                                'no_match',
                                                edgelabel='no_match',
                                                qlabel='accept no match',
                                                color='red',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_glr_glr1_case2_no_match.mask] = 0
            
            M_small_m_nisol_nlr_glr_glr1_case2_match = M_small_m_nisol_nlr_glr_glr1_case2.submask(lofarcat['m_nisol_flag_vc2']==4,
                                                'match',
                                                edgelabel='match',
                                                qlabel='match',
                                                color='blue',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_glr_glr1_case2_match.mask] = 2
            
            M_small_m_nisol_nlr_glr_glr1_case2_deblend = M_small_m_nisol_nlr_glr_glr1_case2.submask(lofarcat['m_nisol_flag_vc2']==5,
                                                'deblend',
                                                edgelabel='deblend',
                                                qlabel='deblend',
                                                color='yellow',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_glr_glr1_case2_deblend.mask] = 4
            
            M_small_m_nisol_nlr_glr_glr1_case2_artefact = M_small_m_nisol_nlr_glr_glr1_case2.submask(lofarcat['m_nisol_flag_vc2']==6,
                                                'artefact',
                                                edgelabel='artefact',
                                                qlabel='artefact',
                                                color='grey',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_glr_glr1_case2_artefact.mask] = -1
            
            
            
            M_small_m_nisol_nlr_glr_glr1_case2_prob = M_small_m_nisol_nlr_glr_glr1_case2.submask(lofarcat['m_nisol_flag_vc2']==0,
                                                'prob',
                                                edgelabel='prob',
                                                qlabel='prob',
                                                #color='grey',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_glr_glr1_case2_prob.mask] = 6
        
        
        
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
        lofarcat[msourceflg][M_small_m_nisol_nlr_glr_glr2_same.mask] = 6
        
        
        # lofarcat['m_nisol_flag_vc2'] = {1,3,4} = lgz, no match, match
        
        
        if step == 3:
            M_small_m_nisol_nlr_glr_glr2_same_lgz = M_small_m_nisol_nlr_glr_glr2_same.submask(lofarcat['m_nisol_flag_vc2'] == 1,
                                                'lgz',
                                                qlabel='lgz',
                                                color='orange',
                                                edgelabel='lgz',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_glr_glr2_same_lgz.mask] = 5
            
            M_small_m_nisol_nlr_glr_glr2_same_no_match = M_small_m_nisol_nlr_glr_glr2_same.submask(lofarcat['m_nisol_flag_vc2'] == 3,
                                                'no_match',
                                                qlabel='accept no match',
                                                color='red',
                                                edgelabel='no_match',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_glr_glr2_same_no_match.mask] = 0
            
            M_small_m_nisol_nlr_glr_glr2_same_match = M_small_m_nisol_nlr_glr_glr2_same.submask(lofarcat['m_nisol_flag_vc2'] == 4,
                                                'match',
                                                qlabel='match',
                                                color='blue',
                                                edgelabel='match',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_glr_glr2_same_match.mask] = 2
            
            
            M_small_m_nisol_nlr_glr_glr2_same_prob = M_small_m_nisol_nlr_glr_glr2_same.submask(lofarcat['m_nisol_flag_vc2'] == 0,
                                                'prob',
                                                qlabel='prob',
                                                #color='blue',
                                                edgelabel='prob',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_glr_glr2_same_prob.mask] = 6
        
        
        M_small_m_nisol_nlr_glr_glr2_diff = M_small_m_nisol_nlr_glr_glr2.submask(lofarcat['Ng_LR_good_unique'] != 1,
                                            'diff',
                                            qlabel='deblend',
                                            color='yellow',
                                            edgelabel='N',
                                            masterlist=masterlist)
        lofarcat[msourceflg][M_small_m_nisol_nlr_glr_glr2_diff.mask] = 4
        
        M_small_m_nisol_nlr_glr_glr3p = M_small_m_nisol_nlr_glr.submask(source_glr3p,
                                            'glr3p',
                                            edgelabel='3+',
                                            qlabel='lgz',
                                            color='green',
                                            masterlist=masterlist)
        lofarcat[msourceflg][M_small_m_nisol_nlr_glr_glr3p.mask] = 5


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
        lofarcat[msourceflg][M_small_m_nisol_nlr_nglr_nglargesep_nnncomplex.mask] = 0
        
        
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
        lofarcat[msourceflg][M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_art.mask] = 0
        
        
        M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart = M_small_m_nisol_nlr_nglr_nglargesep_nncomplex.submask(~nnartefact,
                                            'nnnartefact',
                                            edgelabel='N',
                                            qlabel='visual check',
                                            color=stvcol,
                                            #color='orange',
                                            masterlist=masterlist)
        lofarcat[msourceflg][M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart.mask] = 6
        #1 Complex (lgz)
        #2 complex (lgz - should be done already)
        #3 No match
        #4 artefact
        #5 no match, but NN is artefact
        #6 other (redo)
        #m_nisol_flag_vc1
        if step == 3:
            M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart_lgz = M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart.submask(lofarcat['m_nisol_flag_vc1']==1,
                                                'lgz',
                                                edgelabel='Y',
                                                qlabel='lgz',
                                                color='orange',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart_lgz.mask] = 5 # this was 51
            
            M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart_lgz1 = M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart.submask(lofarcat['m_nisol_flag_vc1']==2,
                                                'lgz1',
                                                edgelabel='Y',
                                                qlabel='lgz (already)',
                                                color='orange',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart_lgz1.mask] = 0  # most of these get associated...
            
            M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart_nlgz = M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart.submask(lofarcat['m_nisol_flag_vc1']>2,
                                                'nlgz',
                                                edgelabel='N',
                                                qlabel='no match',
                                                color='red',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart_nlgz.mask] = 0
            
            
            M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart_prob = M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart.submask(lofarcat['m_nisol_flag_vc1']==0,
                                                'prob',
                                                edgelabel='prob',
                                                qlabel='prob',
                                                #color='red',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_nglr_nglargesep_nncomplex_nart_prob.mask] = 6
        
        M_small_m_nisol_nlr_nglr_glargesep = M_small_m_nisol_nlr_nglr.submask(glarge_sep,
                                            'glargesep',
                                            edgelabel='N',
                                            qlabel='visual check',
                                            color=stvcol,
                                            masterlist=masterlist)
        lofarcat[msourceflg][M_small_m_nisol_nlr_nglr_glargesep.mask] = 6
        
        # lofarcat['m_nisol_flag_vc2'] = {1,3,6} = lgz, no match, artefact
        
        # we've not done visual classification
        if step == 3:
            M_small_m_nisol_nlr_nglr_glargesep_lgz = M_small_m_nisol_nlr_nglr_glargesep.submask(lofarcat['m_nisol_flag_vc2']==1,
                                                'lgz',
                                                edgelabel='lgz',
                                                qlabel='lgz',
                                                color='orange',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_nglr_glargesep_lgz.mask] = 5
            
            
            M_small_m_nisol_nlr_nglr_glargesep_no_match = M_small_m_nisol_nlr_nglr_glargesep.submask(lofarcat['m_nisol_flag_vc2']==3,
                                                'no_match',
                                                edgelabel='no_match',
                                                qlabel='accept no match',
                                                color='red',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_nglr_glargesep_no_match.mask] = 0
            
            
            M_small_m_nisol_nlr_nglr_glargesep_artefact = M_small_m_nisol_nlr_nglr_glargesep.submask(lofarcat['m_nisol_flag_vc2']==6,
                                                'artefact',
                                                edgelabel='artefact',
                                                qlabel='artefact',
                                                color='grey',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_nglr_glargesep_artefact.mask] = -1
            
            M_small_m_nisol_nlr_nglr_glargesep_prob = M_small_m_nisol_nlr_nglr_glargesep.submask(lofarcat['m_nisol_flag_vc2']==0,
                                                'prob',
                                                edgelabel='prob',
                                                qlabel='prob',
                                                #color='grey',
                                                masterlist=masterlist)
            lofarcat[msourceflg][M_small_m_nisol_nlr_nglr_glargesep_prob.mask] = 6
            
            
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
                    #shape='parallelogram'   # end point is a final mask
                    #if t.p < 1.:
                        #label='- {i:n} -\n{lab:s}\n{n:n}\n{p:.2f}%'.format(i=i,lab=t.qlabel,n=t.n,p=t.p)
                    #elif t.p < 10.:
                        #label='- {i:n} -\n{lab:s}\n{n:n}\n{p:.1f}%'.format(i=i,lab=t.qlabel,n=t.n,p=t.p)
                    #else:
                        #label='- {i:n} -\n{lab:s}\n{n:n}\n{p:.0f}%'.format(i=i,lab=t.qlabel,n=t.n,p=t.p)
                        
                    shape='parallelogram'   # end point is a final mask
                    if t.p < 1.:
                        label='- {i:n} -\n{lab:s}\n{n:n}\n{p:.2f}%{mllab:s}'.format(i=i,lab=t.qlabel,n=t.n,p=t.p,mllab=t.mllab)
                    elif t.p < 10.:
                        label='- {i:n} -\n{lab:s}\n{n:n}\n{p:.1f}%{mllab:s}'.format(i=i,lab=t.qlabel,n=t.n,p=t.p,mllab=t.mllab)
                    else:
                        label='- {i:n} -\n{lab:s}\n{n:n}\n{p:.0f}%{mllab:s}'.format(i=i,lab=t.qlabel,n=t.n,p=t.p,mllab=t.mllab)
                        
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
            outname = 'flow_{h}_{v}_step{st:d}_s{s:.0f}_nn{nn:.0f}_msources_{m}_flux{ff:.0f}'.format(h=h,st=step,s=size_large,nn=separation1,v=version,m=mode,ff=fluxcut2)
            A.draw(outname+'.png') # write to file
            A.draw(outname+'.pdf') # write to file
            A.write(outname+'.dot') # write to file

            
        # make a test sample for each final mask
        makesample = 0
        if makesample:
            for t in masterlist:
                if not t.has_children :
                    print(t.name)
                    t.make_sample(lofarcat,Nsample=None)
            
        
        outcat = Table([Column(name='Source_Name', data=lofarcat['Source_Name'][M_small_m_nisol.mask]),
                        Column(name=msourceflg, data=lofarcat[msourceflg][M_small_m_nisol.mask]),
                        Column(name=mcflg, data=lofarcat[mcflg][M_small_m_nisol.mask])])
        #lofarcat.keep_columns(['Source_Name',msourceflg,mcflg])
        #lofarcat[M_small_m_nisol.mask].write(lofar_msource_flowchart_file, overwrite=True)

        #outcat.write(lofar_msource_flowchart_file, overwrite=True)
        
        del lofarcat
        
        update_flags(lofarcat_file_srt, outcat, step)
        
        
        #end of mode
