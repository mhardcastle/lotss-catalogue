import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, join, Column
import astropy.units as u
import astropy.coordinates as ac
import utils.plot_util as pp
import os

from lofar_source_sorter import Mask

if __name__=='__main__':
    
    size_large = 15.           # in arcsec
    separation1 = 45.          # in arcsec
    size_huge = 25.            # in arcsec
    #separation2 = 30.          # in arcsec
    lLR_thresh = 0.36            # LR threshold
    lLR_thresh2 = 0.72            # LR threshold - stricter
    fluxcut = 10               # in mJy
    fluxcut2 = 2.5               # in mJy


    path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'
    lofargcat_file = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.gaus.fixed.fits'
    lofarcat_file = path+'LOFAR_HBA_T1_DR1_catalog_v0.95_masked.srl.fixed.presort.fits'
    psmlcat_file = path+'lofar_pw.fixed.fits'
    psmlgcat_file = path+'lofar_gaus_pw.fixed.fits'

    lofarcat_file_srt = path+'LOFAR_HBA_T1_DR1_catalog_v0.95_masked.srl.fixed.sorted.fits'

    lofar_msource_flowchart_file = path + 'msources/LOFAR_flowchart.fixed.fits'

    # Gaus catalogue
    lofargcat = Table.read(lofargcat_file)
    # only relevant gaussians are in M or C sources
    lofargcat = lofargcat[lofargcat['S_Code'] != 'S']

    # Source catalogue
    lofarcat = Table.read(lofarcat_file_srt)

    # PS ML - matches for sources and gaussians
    psmlcat = Table.read(psmlcat_file)
    psmlgcat = Table.read(psmlgcat_file)

    # from Lara/Philip - classifications of compact isolated m sources
    lofar_msource_flowchart = Table.read(lofar_msource_flowchart_file)


    ### get the panstarrs ML information

    ## join the ps ml cat  - they have identical RA/DEC (source_names were wrong)
    #c = ac.SkyCoord(lofarcat['RA'], lofarcat['DEC'], unit="deg")
    #cpsml = ac.SkyCoord(psmlcat['RA'], psmlcat['DEC'], unit="deg")
    #f_nn_idx,f_nn_sep2d,f_nn_dist3d = ac.match_coordinates_sky(c,cpsml,nthneighbor=1)

    ##psmlcat = psmlcat[f_nn_idx][f_nn_sep2d==0]
    ##lofarcat = lofarcat[f_nn_sep2d==0]

    ## note the large sources are missing from the ML catalogue
    #lrcol = np.zeros(len(lofarcat),dtype=float)
    #lrcol[f_nn_sep2d==0] = psmlcat['lr'][f_nn_idx][f_nn_sep2d==0]

    ##lofarcat.add_column(Column(psmlcat['lr_pc_7th'], 'LR'))
    #lofarcat.add_column(Column(lrcol, 'cLR'))
    #lrcol[np.isnan(lrcol)] = 0
    #lofarcat.add_column(Column(lrcol, 'LR'))
    #lrcol = np.zeros(len(lofarcat),dtype='S19')
    #lrcol[f_nn_sep2d==0] = psmlcat['AllWISE'][f_nn_idx][f_nn_sep2d==0]
    #lofarcat.add_column(Column(lrcol, 'LR_name_wise'))
    #lrcol = np.zeros(len(lofarcat),dtype=int)
    #lrcol[f_nn_sep2d==0] = psmlcat['objID'][f_nn_idx][f_nn_sep2d==0]
    #lofarcat.add_column(Column(lrcol, 'LR_name_ps'))
    #lrcol = np.zeros(len(lofarcat),dtype=float)
    #lrcol[f_nn_sep2d==0] = psmlcat['ra'][f_nn_idx][f_nn_sep2d==0]
    #lofarcat.add_column(Column(lrcol, 'LR_ra'))
    #lrcol = np.zeros(len(lofarcat),dtype=float)
    #lrcol[f_nn_sep2d==0] = psmlcat['dec'][f_nn_idx][f_nn_sep2d==0]
    #lofarcat.add_column(Column(lrcol, 'LR_dec'))


    ## join the ps ml gaus cat  - they have identical RA/DEC (source_names were wrong)
    #cg = ac.SkyCoord(lofargcat['RA'], lofargcat['DEC'], unit="deg")
    #cpsmlg = ac.SkyCoord(psmlgcat['RA'], psmlgcat['DEC'], unit="deg")
    #f_nn_idx_g,f_nn_sep2d_g,f_nn_dist3d_g = ac.match_coordinates_sky(cg,cpsmlg,nthneighbor=1)

    ## note the large sources are missing from the ML catalogue
    #lrgcol = np.zeros(len(lofargcat),dtype=float)
    #lrgcol[f_nn_sep2d_g==0] = psmlgcat['lr'][f_nn_idx_g][f_nn_sep2d_g==0]
    ##lofarcat.add_column(Column(psmlcat['lr_pc_7th'], 'LR'))
    #lofargcat.add_column(Column(lrgcol, 'LR'))
    #lrgcol = np.zeros(len(lofargcat),dtype=float)
    #lrgcol[f_nn_sep2d_g==0] = psmlgcat['ra'][f_nn_idx_g][f_nn_sep2d_g==0]
    #lofargcat.add_column(Column(lrgcol, 'LR_ra'))
    #lrgcol = np.zeros(len(lofargcat),dtype=float)
    #lrgcol[f_nn_sep2d_g==0] = psmlgcat['dec'][f_nn_idx_g][f_nn_sep2d_g==0]
    #lofargcat.add_column(Column(lrgcol, 'LR_dec'))
    #lrgcol = np.zeros(len(lofargcat),dtype='S19')
    #lrgcol[f_nn_sep2d_g==0] = psmlgcat['AllWISE'][f_nn_idx_g][f_nn_sep2d_g==0]
    #lofargcat.add_column(Column(lrgcol, 'LR_name_wise'))
    #lrgcol = np.zeros(len(lofargcat),dtype=int)
    #lrgcol[f_nn_sep2d_g==0] = psmlgcat['objID'][f_nn_idx_g][f_nn_sep2d_g==0]
    #lofargcat.add_column(Column(lrgcol, 'LR_name_ps'))




    #add_G = False   # add the gaussian information
    #lofarcat.add_column(Column(np.ones(len(lofarcat),dtype=int), 'Ng'))
    #lofarcat.add_column(Column(np.ones(len(lofarcat),dtype=float), 'G_LR_max'))
    #lofarcat.add_column(Column(np.ones(len(lofarcat),dtype=int), 'Ng_LR_good'))
    #lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=bool), 'Flag_G_LR_problem'))
    #if add_G:
        #lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=list), 'G_ind'))

    #m_S = lofarcat['S_Code'] =='S'
    #minds = np.where(~m_S)[0]
    #for i,sid in zip(minds, lofarcat['Source_Name'][~m_S]):
        #ig = np.where(lofargcat['Source_Name']==sid)[0]
        #lofarcat['Ng'][i]= len(ig)
        #lofarcat['G_LR_max'][i]= np.nanmax(lofargcat['LR'][ig])
        #igi = np.argmax(lofargcat['LR'][ig])
        ## for now, if one of the gaussian LR is better, take that
        #if lofarcat['G_LR_max'][i] > lofarcat['LR'][i]:
            #lofarcat['LR'][i] = lofarcat['G_LR_max'][i]
            #lofarcat['LR_name_ps'][i] = lofargcat['LR_name_ps'][ig[igi]]
            #lofarcat['LR_name_wise'][i] = lofargcat['LR_name_wise'][ig[igi]]
            #lofarcat['LR_ra'][i] = lofargcat['LR_ra'][ig[igi]]
            #lofarcat['LR_dec'][i] = lofargcat['LR_dec'][ig[igi]]
        ## how many unique acceptable matches are there for the gaussian components
        #matches_ra = np.unique(lofargcat['LR_ra'][ig][np.log10(1+lofargcat['LR'][ig]) > lLR_thresh])
        #n_matches_ra = len(matches_ra)
        #if n_matches_ra > 1:
            #lofarcat['Flag_G_LR_problem'][i] = True
        ## any different to source match
        #if np.sum(matches_ra != lofarcat['LR_ra'][i]):
            #lofarcat['Flag_G_LR_problem'][i] = True
        #lofarcat['Ng_LR_good'][i]= np.nansum(np.log10(1+lofargcat['LR'][ig]) > lLR_thresh)
        
        #if add_G:
            #lofarcat['G_ind'][i]= ig
    #lofarcat['G_LR_max'][m_S] = lofarcat['LR'][m_S]
    #lofarcat['Ng_LR_good'][m_S] = 1*(np.log10(1+lofarcat['LR'][m_S]) > lLR_thresh)



    lofarcat.add_column(Column(-1*np.ones(len(lofarcat),dtype=int), 'M_diag'))

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
    
    m_S_small = (lofarcat['S_Code'] !='S') & (lofarcat['Maj'] <= 15.)
    m_S_small_isol = m_S_small & (lofarcat['NN_sep'] > 45.)
    
    #(lofarcat['Ng_LR_good'][m_S_small_isol] == 0 ) & (source_nlr[m_S_small_isol])
    
    #m_S_small_isol & source_nlr
    
    
    #m_S_small_isol & source_lr
    
    
    masterlist = []
    M_small_m = Mask(m_S_small,
                    'small_m',
                    qlabel='isolated',
                    masterlist=masterlist)
    
    
    
    M_small_m_isol = M_small_m.submask(lofarcat['NN_sep'] <= 45.,
                                        'nisol',
                                        edgelabel='N',
                                        masterlist=masterlist)
    
    M_small_m_nisol = M_small_m.submask(lofarcat['NN_sep'] > 45.,
                                        'isol',
                                        qlabel='LRs >= 0.36??',
                                        edgelabel='Y',
                                        masterlist=masterlist)
    
    M_small_m_isol_lr = M_small_m_isol.submask(source_lr,
                                        'lr',
                                        edgelabel='Y',
                                        qlabel='N(LRg >= 0.36) >= 1?',
                                        masterlist=masterlist)
    
    
    M_small_m_isol_lr_glr = M_small_m_isol_lr.submask(source_glr,
                                        'glr',
                                        edgelabel='Y',
                                        qlabel='N guassian lr?',
                                        masterlist=masterlist)
    
    
    M_small_m_isol_lr_glr_glr1 = M_small_m_isol_lr_glr.submask(source_glr1,
                                        'glr1',
                                        edgelabel='1',
                                        qlabel='same as source?',
                                        masterlist=masterlist)
    
    M_small_m_isol_lr_glr_glr1_ges = M_small_m_isol_lr_glr_glr1.submask(lofarcat['N_G_LR_matchsource'] == 1,
                                        'ges',
                                        edgelabel='Y',
                                        qlabel='accept ML source',
                                        color='blue',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_lr_glr_glr1_ges.mask] = 1
    
    M_small_m_isol_lr_glr_glr1_nges = M_small_m_isol_lr_glr_glr1.submask(lofarcat['N_G_LR_matchsource'] != 1,
                                        'nges',
                                        edgelabel='N',
                                        qlabel='compare LRs and LRg',
                                        masterlist=masterlist)
    
    
    M_small_m_isol_lr_glr_glr1_nges_1 = M_small_m_isol_lr_glr_glr1_nges.submask(lofarcat['G_LR_case1'] == 1,
                                        'case1',
                                        edgelabel='1',
                                        qlabel='accept ML source',
                                        color='blue',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_lr_glr_glr1_nges_1.mask] = 1
    
    
    M_small_m_isol_lr_glr_glr1_nges_2 = M_small_m_isol_lr_glr_glr1_nges.submask(lofarcat['G_LR_case1'] == 2,
                                        'case2',
                                        edgelabel='2',
                                        qlabel='accept ML gaus',
                                        color='blue',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_lr_glr_glr1_nges_2.mask] = 2
    
    
    M_small_m_isol_lr_glr_glr1_nges_3 = M_small_m_isol_lr_glr_glr1_nges.submask(lofarcat['G_LR_case1'] == 3,
                                        'case3',
                                        edgelabel='3',
                                        qlabel='lgz',
                                        color='green',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_lr_glr_glr1_nges_3.mask] = 5
    
    M_small_m_isol_lr_glr_glr3p = M_small_m_isol_lr_glr.submask(source_glr3p,
                                        'glr3p',
                                        qlabel='all G match S?',
                                        edgelabel='3+',
                                        masterlist=masterlist)

    
    M_small_m_isol_lr_glr_glr3p_ges = M_small_m_isol_lr_glr_glr3p.submask(lofarcat['N_G_LR_matchsource'] == lofarcat['Ng_LR_good'],
                                        'ges',
                                        edgelabel='Y',
                                        color='blue',
                                        qlabel='accept ML source',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_lr_glr_glr3p_ges.mask] = 1
    
    M_small_m_isol_lr_glr_glr3p_nges = M_small_m_isol_lr_glr_glr3p.submask(lofarcat['N_G_LR_matchsource'] != lofarcat['Ng_LR_good'],
                                        'nges',
                                        edgelabel='N',
                                        qlabel='deblend',
                                        color='orange',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_lr_glr_glr3p_nges.mask] = 4
    
    
    M_small_m_isol_lr_glr_glr2 = M_small_m_isol_lr_glr.submask(source_glr2,
                                        'compare G/S matches',
                                        edgelabel='2',
                                        masterlist=masterlist)
    
    
    M_small_m_isol_lr_glr_glr2_ges = M_small_m_isol_lr_glr_glr2.submask(lofarcat['N_G_LR_matchsource'] == 2,
                                        'ges',
                                        edgelabel='both G match S',
                                        qlabel='accept ML source',
                                        color='blue',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_lr_glr_glr2_ges.mask] = 1
    
    M_small_m_isol_lr_glr_glr2_nges = M_small_m_isol_lr_glr_glr2.submask(lofarcat['N_G_LR_matchsource'] == 0,
                                        'nges',
                                        edgelabel='neither G match S',
                                        qlabel='lgz',
                                        color='green',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_lr_glr_glr2_nges.mask] = 5
    
    M_small_m_isol_lr_glr_glr2_nges1 = M_small_m_isol_lr_glr_glr2.submask(lofarcat['N_G_LR_matchsource'] == 1,
                                        'nges1',
                                        edgelabel='One G matches S',
                                        qlabel='compare LRsource and LRgaus\'s',
                                        #color='red',
                                        masterlist=masterlist)
    
    M_small_m_isol_lr_glr_glr2_nges1_case1 = M_small_m_isol_lr_glr_glr2_nges1.submask(lofarcat['G_LR_case2'] == 1,
                                        'case1',
                                        edgelabel='1',
                                        qlabel='accept ML source',
                                        color='blue',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_lr_glr_glr2_nges1_case1.mask] = 1
    
    
    M_small_m_isol_lr_glr_glr2_nges1_case2 = M_small_m_isol_lr_glr_glr2_nges1.submask(lofarcat['G_LR_case2'] == 2,
                                        'case2',
                                        edgelabel='2',
                                        qlabel='deblend (direct)',
                                        color='orange',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_lr_glr_glr2_nges1_case2.mask] = 3
    
    
    M_small_m_isol_lr_glr_glr2_nges1_case3 = M_small_m_isol_lr_glr_glr2_nges1.submask(lofarcat['G_LR_case2'] == 3,
                                        'case3',
                                        edgelabel='3',
                                        qlabel='deblend',
                                        color='orange',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_lr_glr_glr2_nges1_case3.mask] = 4
    
    
    M_small_m_isol_lr_nglr = M_small_m_isol_lr.submask(source_nglr,
                                        'nglr',
                                        edgelabel='N',
                                        qlabel='LRs >= 10*0.36?',
                                        masterlist=masterlist)
    
    M_small_m_isol_lr_nglr_lr2 = M_small_m_isol_lr_nglr.submask(source_lr2,
                                        'slr',
                                        edgelabel='Y',
                                        qlabel='accept ML source',
                                        color='blue',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_lr_nglr_lr2.mask] = 1
    
    M_small_m_isol_lr_nglr_nlr2 = M_small_m_isol_lr_nglr.submask(~source_lr2,
                                        'nslr',
                                        edgelabel='N',
                                        qlabel='lgz',
                                        color='green',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_lr_nglr_nlr2.mask] = 5
    
    M_small_m_isol_nlr = M_small_m_isol.submask(source_nlr,
                                        'nlr',
                                        edgelabel='N',
                                        qlabel='N(LRg >= 0.36) >= 1?',
                                        masterlist=masterlist)
    

    M_small_m_isol_nlr_glr = M_small_m_isol_nlr.submask(source_glr,
                                        'glr',
                                        edgelabel='Y',
                                        qlabel='N(LRg >= 0.36)?',
                                        masterlist=masterlist)
    
    
    M_small_m_isol_nlr_glr_glr1 = M_small_m_isol_nlr_glr.submask(source_glr1,
                                        'glr1',
                                        edgelabel='1',
                                        qlabel='LRgaus > 10*0.36 & \nMaj_gaus< 10"?',
                                        masterlist=masterlist)
    
    M_small_m_isol_nlr_glr_glr1_case1 = M_small_m_isol_nlr_glr_glr1.submask(lofarcat['G_LR_case3']==1,
                                        'case1',
                                        edgelabel='Y',
                                        qlabel='accept ML gaus',
                                        color='blue',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_nlr_glr_glr1_case1.mask] = 2
    
    M_small_m_isol_nlr_glr_glr1_case2 = M_small_m_isol_nlr_glr_glr1.submask(lofarcat['G_LR_case3']==2,
                                        'case2',
                                        edgelabel='N',
                                        qlabel='lgz',
                                        color='green',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_nlr_glr_glr1_case2.mask] = 5
    
    M_small_m_isol_nlr_glr_glr2 = M_small_m_isol_nlr_glr.submask(source_glr2,
                                        'glr2',
                                        qlabel='match same gal?',
                                        edgelabel='2',
                                        masterlist=masterlist)
    
    
    M_small_m_isol_nlr_glr_glr2_same = M_small_m_isol_nlr_glr_glr2.submask(lofarcat['Ng_LR_good_unique'] == 1,
                                        'same',
                                        qlabel='lgz',
                                        color='green',
                                        edgelabel='Y',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_nlr_glr_glr2_same.mask] = 5
    
    
    M_small_m_isol_nlr_glr_glr2_diff = M_small_m_isol_nlr_glr_glr2.submask(lofarcat['Ng_LR_good_unique'] != 1,
                                        'diff',
                                        qlabel='deblend',
                                        color='orange',
                                        edgelabel='N',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_nlr_glr_glr2_diff.mask] = 4
    
    M_small_m_isol_nlr_glr_glr3p = M_small_m_isol_nlr_glr.submask(source_glr3p,
                                        'glr3p',
                                        edgelabel='3+',
                                        qlabel='lgz',
                                        color='green',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_nlr_glr_glr3p.mask] = 5


    M_small_m_isol_nlr_nglr = M_small_m_isol_nlr.submask(source_nglr,
                                        'nglr',
                                        edgelabel='N',
                                        qlabel='max(sep_g) < 15"?',
                                        masterlist=masterlist)
    
    
    
    M_small_m_isol_nlr_nglr_nglargesep = M_small_m_isol_nlr_nglr.submask(~glarge_sep,
                                        'nglargesep',
                                        edgelabel='Y',
                                        qlabel='NNsep<100" & NNsize>10"?',
                                        masterlist=masterlist)
    
    
    M_small_m_isol_nlr_nglr_nglargesep_nnncomplex = M_small_m_isol_nlr_nglr_nglargesep.submask(~nncomplex,
                                        'nnncomplex',
                                        edgelabel='N',
                                        qlabel='accept no match',
                                        color='red',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_nlr_nglr_nglargesep_nnncomplex.mask] = 0
    
    
    M_small_m_isol_nlr_nglr_nglargesep_nncomplex = M_small_m_isol_nlr_nglr_nglargesep.submask(nncomplex,
                                        'nncomplex',
                                        edgelabel='Y',
                                        qlabel='lgz',
                                        color='green',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_nlr_nglr_nglargesep_nncomplex.mask] = 5
    
    
    
    M_small_m_isol_nlr_nglr_glargesep = M_small_m_isol_nlr_nglr.submask(glarge_sep,
                                        'glargesep',
                                        edgelabel='N',
                                        qlabel='lgz',
                                        color='green',
                                        masterlist=masterlist)
    lofarcat['M_diag'][M_small_m_isol_nlr_nglr_glargesep.mask] = 5
    
    
    # make flowchart from list of masks
    plot_flowchart = True
    plot_verbose = False
    try:
        import pygraphviz as pgv
    except ImportError:
        print 'no pygraphviz; cannot make visual flowchart'
        plot_flowchart = False
    if plot_flowchart:

        PW = 1500.

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
            print(A.string()) # print dot file to standard output

        # make the flowchart
        #Optional prog=['neato'|'dot'|'twopi'|'circo'|'fdp'|'nop']
        #neato, dot, twopi, circo, fdp, nop, wc, acyclic, gvpr, gvcolor, ccomps, sccmap, tred, sfdp.
        A.layout('dot') # layout with dot
        A.draw('flow_s{s:.0f}_nn{nn:.0f}_msources_nisol.png'.format(s=size_large,nn=separation1)) # write to file
        A.write('flow_s{s:.0f}_nn{nn:.0f}_msources_nisol.dot'.format(s=size_large,nn=separation1)) # write to file

        