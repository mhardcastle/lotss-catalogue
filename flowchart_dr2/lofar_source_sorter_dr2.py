'''
lofar_source_sorter.py
the main code for the decision tree
'''
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from astropy.table import Table, join, Column
import astropy.units as u
import astropy.coordinates as ac

import utils.plotting as pp


class Mask:
    '''Mask store a boolean mask and associated information necessary for the flowchart
    mask - boolean mask to apply
    trait - short description
    label - long detailed description (default: same as trait)
    level - level of nested masks (defatul: 0)
    verbose - print some of output (default: True)
    masterlist - add output mask to list (default: None)
    qlabel - for flowchart: what question gets asked of this level (default: None)
    color - for flowchart: color to plot the mask (default: None=black)
    '''
    
    def __init__(self, mask, trait, label=None, level=0, verbose=True, masterlist=None, qlabel=None, color=None, maskflag=None, maskflagstr=''):
        self.mask = mask
        if qlabel is not None :
            self.qlabel = qlabel
        else:
            self.qlabel = label
        if label is None:
            label = trait
        self.label = label
        if isinstance(trait,str):
            self.traits = list([trait])
        else:
            self.traits = list(trait)
            
        self.name = '_'.join(self.traits)
            
        self.color = color
        self.level = level
        
        self.N = self.total()
        self.n = self.msum()
        self.f = self.fraction()
        self.p = self.percent()
        
        self.has_children = False
        self.has_parent = False
        
        self.maskflagstr = maskflagstr
        self.maskflag = maskflag
        if maskflag is not None:
            self.mask2 = maskflag  #[self.mask]
            self.set_mllab()
        else:
            self.mask2 = self.mask
            self.set_mllab()
        
        self.Nchildren = 0
        self.children = None
        self.parent = None
        
        if masterlist is not None:
            masterlist.append(self)
        
        if verbose:
            self.print_frac()
        
        return
    
    def set_mllab(self):
        if self.mask2 is None:
            self.mllab = ''
        else:
            n1 = np.sum(self.mask & (self.mask2==1))
            n2 = np.sum(self.mask & (self.mask2==0))
            Ni = np.sum(self.mask)
            p1 = 100.*n1/Ni
            p2 = 100.*n2/Ni
            mlstr = self.maskflagstr
            self.mllab='\n{mlstr:s} {n1:d}:{n2:d}  ({p1:.0f}:{p2:.0f}%)\n'.format(mlstr=mlstr,n1=n1,n2=n2,p1=p1,p2=p2)
    
    def percent(self):
        return 100.*np.sum(self.mask)/self.N
    
    def fraction(self):
        return 1.*np.sum(self.mask)/self.N

    def msum(self):
        return np.sum(self.mask)
    
    def total(self):
        return len(self.mask)
    
    def print_frac(self, vformat=True):
        '''vformat = True will print with formatted spaces indicative of the hierarchical structure
        '''
        if vformat and self.level > 0:
            vv = ' '*self.level + '-'*self.level
        else:
            vv = ' '
        print('{n:6d} ({f:5.1f}%){vv:s}{label:s}'.format(vv=vv, n=self.n, f=self.p, label=self.label))
        
    def __str__(self):
        return self.name
        
    def submask(self, joinmask, newtrait, label=None, edgelabel='Y', verbose=True, qlabel=None, masterlist=None, color=None):
        '''create a new submask based on this instance -- join masks with AND
        # qlabel  is the question that will be asked
        # edgelabel is the answer to the question asked to get here
        '''
        newmask = self.mask & joinmask
        newmaskflag = self.maskflag 
        newtraits = list(self.traits)  # copy list of traits - lists are mutable!!
        newtraits.append(newtrait)     # append new trait onto copy
        newlevel = self.level + 1
        
        
        if label is None:
            label = newtrait
        
        childmask = Mask(newmask, newtraits, label, level=newlevel, masterlist=masterlist, verbose=verbose, qlabel=qlabel, color=color, maskflag=newmaskflag, maskflagstr=self.maskflagstr)  
        
        childmask.has_parent = True
        childmask.parent = self
        
        childmask.edgelabel = edgelabel  
        
        if not self.has_children:
            self.has_children = True
            self.children = [childmask]
            self.Nchildren = 1
        else:
            newchildren = list(self.children)  # copy list of traits - lists are mutable!!
            newchildren.append(childmask)
            self.children = newchildren
            self.Nchildren = len(newchildren)
            
        return childmask
    
    # make sample files
    def make_sample(self, cat, Nsample=250):
        '''create a random subsample of the masked catalogue 'cat'
        '''
        
        t = cat[self.mask]
        if Nsample is None:
            Nsample = len(t)
        Nsample = np.min((Nsample, len(t)))
        if Nsample ==0 : return
        if Nsample < len(t):
            t = t[np.random.choice(np.arange(len(t)), Nsample, replace=False)]
        fitsname = 'sample_'+self.name+'.fits'
        if os.path.exists(fitsname):
            os.remove(fitsname)
        if not os.path.exists('flowsamples'):
            os.mkdir('flowsamples')
        t.write('flowsamples/'+fitsname, overwrite=True)
        
        return
    
    def is_disjoint(self, othermask):
        
        assert isinstance(othermask, Mask), 'need to compare to another Mask instance'
        
        if np.sum((self.mask) & (othermask.mask)) == 0:
            return True
        else:
            return False
        
        return
        
def Masks_disjoint_complete(masklist):
    '''test whether a list of masks is disjoint and complete
    '''
    Z = np.zeros(len(masklist[0].mask), dtype=bool)
    O = np.ones(len(masklist[0].mask), dtype=bool)
    for t in masklist:
        Z = Z & t.mask
        O = O | t.mask
    
    return np.all(O) and np.all(~Z)


if __name__=='__main__':


    step = int(sys.argv[1])
    if step not in  [1,2,3,4]:
        print('unknown step',step)
        sys.exit(1)
    # step 1 is first assign FC_flag1
    # step 2 is after running msource stuff
    # step 3 is after running msource stuff and inputs from ML
    # step 4 is after visual classification if we do that this time
    
        
    size_large = 15.           # in arcsec
    separation1 = 45.          # in arcsec
    size_huge = 25.            # in arcsec
    #separation2 = 30.          # in arcsec
    lLR_thresh = 0.404            # LR threshold
    fluxcut = 8.               # in mJy
    fluxcut2 = 2.5               # in mJy

    ### Required INPUTS
    # lofar source catalogue, gaussian catalogue and ML catalogues for each

    path = '/data2/wwilliams/projects/lofar_surveys/LoTSS-DR2-Feb2020/'
    
    if step == 1:
        # lr infor is in the catalogue - this makes things easier
        lofarcat_file = path+'LoTSS_DR2_rolling.srl_0h.lr.presort.fits'
        lofargcat_file = path+'LoTSS_DR2_rolling.gaus_0h.lr.fits'
    elif step >= 2:
        lofarcat_file = path+'LoTSS_DR2_rolling.srl_0h.sorted_step1.fits'
        lofargcat_file = path+'LoTSS_DR2_rolling.gaus_0h.lr.fits'
        


    # output file - with updated columns
    lofarcat_file_srt = path+'LoTSS_DR2_rolling.srl_0h.sorted_step{st}.fits'.format(st=step)



    # Gaus catalogue
    lofargcat = Table.read(lofargcat_file)
    # only relevant gaussians are in M or C sources
    lofargcat = lofargcat[lofargcat['S_Code'] != 'S']

    # Source catalogue
    lofarcat = Table.read(lofarcat_file)

    # these two added with add_gaus_info.py
    needgaus = False
    for gaus_col in ['Ng', 'G_max_sep', 'G_LR_max', 'Ng_LR_good','Ng_LR_good_unique','N_G_LR_matchsource','Flag_G_LR_problem']:
        if gaus_col not in lofarcat.colnames:
            print('missing gaussian column: ',gaus_col)
            needgaus = True
    if needgaus:
        raise  RuntimeError('need the Gaussian count and max separation and LR information -- run add_gaus_info.py')
    
    
       
    
    if step == 1:
        print('msources have not yet been handled, but we need to give them FC_flags')
    elif step == 2:
        print('msources have been handled, now to give them the right ID flags')
        if 'msource1_flag' not in lofarcat.colnames:
            raise  RuntimeError('need the msource1_flag -- run get_msource_flag.py')
        if 'msource2_flag' not in lofarcat.colnames:
            raise  RuntimeError('need the msource2_flag -- run get_msource_flag.py')
        fixlist = ['ILTJ111611.15+493234.8', 'ILTJ121557.43+512418.5', 'ILTJ121956.33+473756.2', 'ILTJ124047.70+500956.4', 'ILTJ131951.11+531748.0', 'ILTJ135441.11+541646.4']
        for s in fixlist:
            lofarcat['msource1_flag'][lofarcat['Source_Name']==s] = 1
        fixlist = ['ILTJ140430.41+560657.9'] 
        for s in fixlist:
            lofarcat['msource1_flag'][lofarcat['Source_Name']==s] = 2
    elif step == 3:
        print('need ML outputs')
        print('machine learning has been run, now to give them the right ID flags')
        if 'ML_flag' not in lofarcat.colnames:
            raise  RuntimeError('need the ML_flag -- run get_ml_flags.py')  ## 1 for LR, 0 for LGZ
    else:
        print('step {s} not defined, quitting'.format(s=step))
        sys.exit(1)



    # this is easy to run...
    # get the large 2masx sources (must run match_2masx for these)
    if '2MASX_match_large' not in lofarcat.colnames:
        raise  RuntimeError('need the 2masx information -- run match_2masx.py')
    big2masx = lofarcat['2MASX_match_large']

    # step 1 happens before we have visual classifications - if any ever for dr2
    if step in [1,2,3]:
        # we start with no flagging - visual, other
        cleanflag = np.ones(len(lofarcat),dtype=bool)
        cleanflag0 = np.zeros(len(lofarcat),dtype=bool)
        clustered_flag = cleanflag0
        Lclustered_flag = cleanflag0
        huge_faint_flag = cleanflag0
        nhuge_faint_flag = cleanflag0
        nhuge_2masx_flag = cleanflag0
        double_flag = cleanflag0
        double_flag2 = cleanflag0
        mnisol_flag2 = cleanflag0
        mnisol_flag1 = cleanflag0
        #big2masx = cleanflag0
        artefact_flag = cleanflag0
        edge_flag = cleanflag0
        # combine the artefact flags
        # artefacts have been identified through various routes of visual checking
        Artefact_flag = (artefact_flag == 1)
        
        if '2MASX' not in lofarcat.colnames:
            lofarcat.add_column(Column(data=cleanflag0, name='2MASX'))
        if 'double_flag' not in lofarcat.colnames:
            lofarcat.add_column(Column(data=cleanflag0, name='double_flag'))
        
    elif step == 4:
        # get the visual flags (must run get_visual_flags for these after doing visual confirmation - classify_*.py)
        if 'clustered_flag' not in lofarcat.colnames:
            raise  RuntimeError('need the visual flag information for the clustered sources')
        clustered_flag = lofarcat['clustered_flag']
            
        if 'Lclustered_flag' not in lofarcat.colnames:
            raise  RuntimeError('need the visual flag information for the large faint clustered sources')
        Lclustered_flag = lofarcat['Lclustered_flag']
            
        if 'huge_faint_flag' not in lofarcat.colnames:
            raise  RuntimeError('need the visual flag information for the huge faint sources')
        huge_faint_flag = lofarcat['huge_faint_flag']

        if 'nhuge_faint_flag' not in lofarcat.colnames:
            raise  RuntimeError('need the visual flag information for the huge faint sources')
        nhuge_faint_flag = lofarcat['nhuge_faint_flag']
            
        if 'nhuge_2masx_flag' not in lofarcat.colnames:
            raise  RuntimeError('need the visual flag information for the large_nhuge_2masx sources')
        nhuge_2masx_flag = lofarcat['nhuge_2masx_flag']

        if 'double_flag' not in lofarcat.colnames:
            raise  RuntimeError('need the visual flag information for the double sources')
        double_flag = lofarcat['double_flag']
                    
        if 'double_flag2' not in lofarcat.colnames:
            raise  RuntimeError('need the visual flag information for the double sources')
        double_flag2 = lofarcat['double_flag2']
                    
        if 'm_nisol_flag_vc2' not in lofarcat.colnames:
            raise  RuntimeError('need the  visual flag information for the m non-isol sources')
        mnisol_flag2 = lofarcat['m_nisol_flag_vc2']
                    
        if 'm_nisol_flag_vc1' not in lofarcat.colnames:
            raise  RuntimeError('need the  visual flag information for the m non-isol sources')
        mnisol_flag1 = lofarcat['m_nisol_flag_vc1']


        ## get artefact information (must run find_artefacts for these)
        if 'artefact_flag' not in lofarcat.colnames:
            raise  RuntimeError('need the artefact information')
        artefact_flag = lofarcat['artefact_flag']
            
        ## get edge_flag information (must run flag_edge_sources.py and get_visual_flags for these)
        if 'edge_flag' not in lofarcat.colnames:
            raise  RuntimeError('need the edge flags from get_visual_flags')
        edge_flag = lofarcat['edge_flag']
            
            
        # combine the artefact flags
        # artefacts have been identified through various routes of visual checking
        Artefact_flag = (artefact_flag == 1) | (huge_faint_flag ==4) | (nhuge_2masx_flag==4) | (Lclustered_flag == 1) | (clustered_flag == 1) | (nhuge_faint_flag==5) | (edge_flag==True) | (mnisol_flag2 == 6) | (double_flag == 3) | (double_flag2 == 3)


    if 'Artefact_flag' not in lofarcat.colnames:
        lofarcat.add_column(Column(Artefact_flag, 'Artefact_flag'))
    else:
        lofarcat['Artefact_flag'] = Artefact_flag

    if 'ID_flag' not in lofarcat.colnames:
        lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=int),'ID_flag'))
    if step in [1,2]:
        lofarcat['ID_flag'] = -1*np.ones(len(lofarcat),dtype=int)
    
    if 'LGZ_flag' not in lofarcat.colnames:
        lofarcat.add_column(Column(np.zeros(len(lofarcat),dtype=int),'LGZ_flag'))
    if step in [1,2]:
        lofarcat['LGZ_flag'] = np.zeros(len(lofarcat),dtype=int)
        
    

    #############################################################################


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
        
    # now exclude artefacts - just put them far away always at the south pole
    dec = lofarcat['DEC'].copy()
    dec[Artefact_flag] = -90
    cclean = ac.SkyCoord(lofarcat['RA'], dec, unit="deg")
    f_nnc_idx,f_nnc_sep2d,_ = ac.match_coordinates_sky(cclean,cclean,nthneighbor=2)
    f_nnc4_idx,f_nnc4_sep2d,_ = ac.match_coordinates_sky(cclean,cclean,nthneighbor=4)
    f_nnc5_idx,f_nnc5_sep2d,_ = ac.match_coordinates_sky(cclean,cclean,nthneighbor=5)

    #'clean' nearest neighbour
    if 'NNC_sep' not in lofarcat.colnames:
        lofarcat.add_column(Column(lofarcat['LR'][f_nnc_idx], 'NNC_LR'))
        lofarcat.add_column(Column(f_nnc_sep2d.to(u.arcsec).value, 'NNC_sep'))
        lofarcat.add_column(Column(f_nnc_idx, 'NNC_idx'))
        lofarcat.add_column(Column(f_nnc5_sep2d.to(u.arcsec).value, 'NNC5_sep'))
        lofarcat.add_column(Column(f_nnc4_sep2d.to(u.arcsec).value, 'NNC4_sep'))
        lofarcat.add_column(Column(lofarcat['Total_flux'][f_nnc_idx], 'NNC_Total_flux'))
        lofarcat.add_column(Column(lofarcat['Total_flux']/lofarcat['NNC_Total_flux'], 'NNC_Frat'))
        lofarcat.add_column(Column(lofarcat['Maj'][f_nnc_idx], 'NNC_Maj'))





    ########################################################

    Ncat = len(lofarcat)

    #m_all = lofarcat['RA'] > -1

    masterlist = []

    M_all = Mask(lofarcat['RA'] > -1,
                    'all',
                    qlabel='artefact?\n(visual confirmation)',
                    masterlist=masterlist,
                    maskflag=lofarcat['ML_flag'], maskflagstr='LR:LGZ')

    # artefacts
    M_all_artefact = M_all.submask(artefact_flag,
                        'artefact',
                        'Artefact\n(visually confirmed)',
                        edgelabel='Y',
                        color='gray',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_all_artefact.mask] = -1

    # sources 
    M_all_clean = M_all.submask(~artefact_flag,
                        'src',
                        'Clean',
                        edgelabel='N',
                        qlabel='Huge 2MASX?\n(r(2MASX)>60")',
                        masterlist=masterlist)


    # big optical gal 
    M_all_biggal = M_all_clean.submask((big2masx),
                        'big2MASX',
                        'Huge 2MASX source\n(64 r(2MASX)>60" galaxies)',
                        edgelabel='Y',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_all_biggal.mask] = 2

    # sources - clean
    M_all_clean2 = M_all_clean.submask(~big2masx,
                        'clean',
                        'Clean2',
                        edgelabel='N',
                        qlabel='Large?\n(s>{s:.0f}")'.format(s=size_large),
                        masterlist=masterlist)


    # large 
    M_large = M_all_clean2.submask(lofarcat['Maj'] > size_large,
                        'large',
                        'large (s>{s:.0f}")'.format(s=size_large),
                        qlabel='Bright?\n(S>{f:.0f} mJy)'.format(f=fluxcut, s=size_large),
                        masterlist=masterlist)

    # large bright
    M_large_bright = M_large.submask(lofarcat['Total_flux'] > fluxcut,
                        'bright',
                        'large (s>{s:.0f}") & bright (S>{f:.0f} mJy)'.format(f=fluxcut, s=size_large),
                        qlabel='LGZ',
                        color='green',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_large_bright.mask] = 31
    lofarcat['LGZ_flag'][M_large_bright.mask] = 1

    # large faint
    M_large_faint = M_large.submask(lofarcat['Total_flux'] <= fluxcut,
                        'faint',
                        'large (s>{s:.0f}") & faint (S<={f:.0f} mJy)'.format(f=fluxcut, s=size_large),
                        edgelabel='N',
                        qlabel='Visual sorting',
                        color='orange',
                        masterlist=masterlist)
    
    lofarcat['ID_flag'][M_large_faint.mask] = 5  #5 is tbd

    # huge faint 
    huge = (lofarcat['Maj'] > 2*size_large)
    
    # artefacts
    lf_artefacts =  (huge & (huge_faint_flag ==4)) | \
        (lofarcat['2MASX'] & (nhuge_2masx_flag==4)) | \
        ((lofarcat['NN4_sep'] <= separation1+size_large) & (Lclustered_flag == 1)) | \
        (nhuge_faint_flag == 5)
    
    # complex
    lf_complex = (huge & (huge_faint_flag ==1)) | \
        (lofarcat['2MASX'] & (nhuge_2masx_flag==2)) | \
        ((lofarcat['NN4_sep'] <= separation1+size_large) & (Lclustered_flag == 2)) | \
        (nhuge_faint_flag == 1)
    
    # complex-zoom
    lf_complex_zoom = (nhuge_faint_flag == 4)
    
    # bright galaxy
    lf_bright =  (lofarcat['2MASX'] & (nhuge_2masx_flag==1)) | \
        (huge & (huge_faint_flag ==2))
    
    # no match possible
    lf_nomatch =  ( huge & (huge_faint_flag ==3)) | \
        (nhuge_faint_flag == 3)
    
    # good lr
    lf_match =   (nhuge_faint_flag == 2)
    
    ## we have outputs from machine learning ... can go here
    #if step == 3:
        #M_large_faint_lgz = M_large_faint.submask(lf_artefacts,
                            #'artefact',
                            #edgelabel='N(r)',
                            #qlabel='artefact\n(visually confirmed)',
                            #color='gray',
                            #masterlist=masterlist)
        #lofarcat['ID_flag'][M_large_faint_artefact.mask] = -1
        
    # if we have outputs from visual sorting we can go here
    if step == 4:
        M_large_faint_artefact = M_large_faint.submask(lf_artefacts,
                            'artefact',
                            edgelabel='N(r)',
                            qlabel='artefact\n(visually confirmed)',
                            color='gray',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_large_faint_artefact.mask] = -1
        
        M_large_faint_complex = M_large_faint.submask(lf_complex,
                            'complex',
                            edgelabel='Y(*)',
                            qlabel='complex\n(LGZ)',
                            color='green',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_large_faint_complex.mask] = 3
        lofarcat['LGZ_flag'][M_large_faint_complex.mask] = 2
        
        M_large_faint_complex_zoom = M_large_faint.submask(lf_complex_zoom,
                            'lgzz',
                            edgelabel='complex-zoom',
                            qlabel='LGZ-zoom',
                            color='green',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_large_faint_complex_zoom.mask] = 3
        lofarcat['LGZ_flag'][M_large_faint_complex_zoom.mask] = 20
        
        
        M_large_faint_match = M_large_faint.submask(lf_match,
                            'match',
                            edgelabel='Y(m)',
                            qlabel='Accept ML\n(visually confirmed)',
                            color='blue',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_large_faint_match.mask] = 1
        
        M_large_faint_nomatch = M_large_faint.submask(lf_nomatch,
                            'nomatch',
                            edgelabel='Y(nm)',
                            qlabel='No match possible\n(visually confirmed)',
                            color='red',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_large_faint_nomatch.mask] = 4
            
        M_large_faint_2masx = M_large_faint.submask(lf_bright,
                            '2masx',
                            edgelabel='Y',
                            qlabel='bright galaxy\n(visually confirmed)',
                            color='blue',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_large_faint_2masx.mask] = 2
        

    #####

    # compact 
    M_small = M_all_clean2.submask(lofarcat['Maj'] <= size_large,
                        'small',
                        'compact (s<{s:.0f}")'.format(s=size_large),
                        edgelabel='N',
                        qlabel='Isolated?\n(NN>{nn:.0f}")'.format(nn=separation1),
                        masterlist=masterlist)
    # compact isolated
    M_small_isol = M_small.submask(lofarcat['NN_sep'] > separation1,
                        'isol',
                        'compact isolated (s<{s:.0f}", NN>{nn:.0f}")'.format(s=size_large, nn=separation1),
                        qlabel='S?',
                        masterlist=masterlist)


    # compact isolated
    M_small_isol_S = M_small_isol.submask(lofarcat['S_Code'] == 'S',
                        'S',
                        'compact isolated (s<{s:.0f}", NN>{nn:.0f}") S'.format(s=size_large, nn=separation1),
                        qlabel='LR > {l:.2f}?'.format(l=lLR_thresh),
                        masterlist=masterlist)


    # compact isolated good lr
    M_small_isol_S_lr = M_small_isol_S.submask(lofarcat['LR'] >= lLR_thresh,
                        'lr',
                        'compact isolated good LR (s<{s:.0f}", NN>{nn:.0f}")'.format(s=size_large, nn=separation1),
                        color='blue',
                        qlabel='Accept LR',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_isol_S_lr.mask] = 1


    # compact isolated badd lr
    M_small_isol_S_nlr = M_small_isol_S.submask(lofarcat['LR'] < lLR_thresh,
                        'nlr',
                        'compact isolated bad LR (s<{s:.0f}", NN>{nn:.0f}")'.format(s=size_large, nn=separation1),
                        edgelabel='N',
                        color='red',
                        qlabel='Accept no LR',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_isol_S_nlr.mask] = 1


    # compact isolated nS
    M_small_isol_nS = M_small_isol.submask(lofarcat['S_Code'] != 'S',
                        'nS',
                        'compact isolated (s<{s:.0f}", NN>{nn:.0f}") !S'.format(s=size_large, nn=separation1),
                        edgelabel='N',
                        color='gray',
                        qlabel='msource',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_isol_nS.mask] = 5  # - these should all be overwritten below, but leave it here to doublecheck (there should be no 5's in the end)
    if step == 2:
        # the ID_flags here depend on the output of Philips flow chart:
        lofarcat['ID_flag'][M_small_isol_nS.mask & (lofarcat['msource1_flag'] == 0) ] = 12
        lofarcat['ID_flag'][M_small_isol_nS.mask & (lofarcat['msource1_flag'] == 1) ] = 12
        lofarcat['ID_flag'][M_small_isol_nS.mask & (lofarcat['msource1_flag'] == 2) ] = 12
        lofarcat['ID_flag'][M_small_isol_nS.mask & (lofarcat['msource1_flag'] == 3) ] = 41
        lofarcat['ID_flag'][M_small_isol_nS.mask & (lofarcat['msource1_flag'] == 4) ] = 42
        lofarcat['ID_flag'][M_small_isol_nS.mask & (lofarcat['msource1_flag'] >= 5) ] = 32
        lofarcat['LGZ_flag'][M_small_isol_nS.mask & (lofarcat['msource1_flag'] >= 5) ] = 20

    



    # compact not isolated
    M_small_nisol = M_small.submask(lofarcat['NN_sep'] <= separation1,
                        'nisol',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}")'.format(s=size_large, nn=separation1),
                        edgelabel='N',
                        qlabel='Clustered?\n(NN5<{nn:.0f}")'.format(s=size_large, nn=separation1),
                        masterlist=masterlist)
    
    if step == 4:
        M_small_nisol_artefact = M_small_nisol.submask((lofarcat['NN5_sep'] <= separation1) & (clustered_flag == 1),
                            'artefact',
                            edgelabel='N(r)',
                            qlabel='artefact\n(visually confirmated)',
                            color='gray',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_small_nisol_artefact.mask] = -1
        
        M_small_nisol_complex = M_small_nisol.submask((lofarcat['NN5_sep'] <= separation1) & (clustered_flag == 2),
                            'complex',
                            edgelabel='Y',
                            qlabel='complex\n(LGZ)',
                            color='green',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_small_nisol_complex.mask] = 5
        lofarcat['LGZ_flag'][M_small_nisol_complex.mask] = 2
        
        M_small_nisol_nclustered = M_small_nisol.submask((lofarcat['NN5_sep'] > separation1) | ((lofarcat['NN5_sep'] <= separation1) & (clustered_flag == 3)),
                            'nclustered',
                            'compact not isolated (s<{s:.0f}", NN5>{nn:.0f}")'.format(s=size_large, nn=separation1),
                            edgelabel='N',
                            qlabel='S?',
                            masterlist=masterlist)
        
    elif step in [1,2,3]:
        # we have no clustered flag
        M_small_nisol_tbc = M_small_nisol.submask((lofarcat['NN5_sep'] <= separation1),
                            'tbc',
                            edgelabel='Y',
                            qlabel='clustered\n(Visual sorting)',
                            color='orange',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_small_nisol_tbc.mask] = 5

        M_small_nisol_nclustered = M_small_nisol.submask((lofarcat['NN5_sep'] > separation1),
                            'nclustered',
                            'compact not isolated (s<{s:.0f}", NN5>{nn:.0f}")'.format(s=size_large, nn=separation1),
                            edgelabel='N',
                            qlabel='S?',
                            masterlist=masterlist)

    # compact not isolated, nnlarge
    M_small_nisol_nclustered_nS = M_small_nisol_nclustered.submask(lofarcat['S_Code'] != 'S',
                        'nS',
                        edgelabel='N',
                        qlabel='msource'.format(l=lLR_thresh),
                        color='gray',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_nclustered_nS.mask] = 5  # - these should all be overwritten below, but leave it here to doublecheck (there should be no 5's in the end)

    if step == 2:
        # the ID_flags here depend on the output of Philips flow chart:
        lofarcat['ID_flag'][M_small_nisol_nclustered_nS.mask & (lofarcat['msource2_flag'] == 0) ] = 12
        lofarcat['ID_flag'][M_small_nisol_nclustered_nS.mask & (lofarcat['msource2_flag'] == 1) ] = 12
        lofarcat['ID_flag'][M_small_nisol_nclustered_nS.mask & (lofarcat['msource2_flag'] == 2) ] = 12
        lofarcat['ID_flag'][M_small_nisol_nclustered_nS.mask & (lofarcat['msource2_flag'] == 3) ] = 41
        lofarcat['ID_flag'][M_small_nisol_nclustered_nS.mask & (lofarcat['msource2_flag'] == 4) ] = 42
        lofarcat['ID_flag'][M_small_nisol_nclustered_nS.mask & (lofarcat['msource2_flag'] >= 5) ] = 32
        lofarcat['LGZ_flag'][M_small_nisol_nclustered_nS.mask & (lofarcat['msource2_flag'] >= 5) ] = 20

    

    # compact not isolated, nnlarge
    M_small_nisol_nclustered_S = M_small_nisol_nclustered.submask(lofarcat['S_Code'] == 'S',
                        'S',
                        edgelabel='Y',
                        qlabel='LR > {l:.2f}?'.format(l=lLR_thresh),
                        masterlist=masterlist)

    # compact not isolated, nnsmall, lr
    M_small_nisol_nclustered_S_lr = M_small_nisol_nclustered_S.submask(lofarcat['LR'] >= lLR_thresh,
                        'lr',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), good LR'.format(s=size_large, nn=separation1),
                        edgelabel='Y',
                        qlabel='NN LR > {l:.2f}?'.format(l=lLR_thresh),
                        masterlist=masterlist)

    # compact not isolated, nnsmall, lr, NNlr
    M_small_nisol_nclustered_S_lr_NNlr = M_small_nisol_nclustered_S_lr.submask(lofarcat['NN_LR'] >= lLR_thresh,
                        'NNlr',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), good LR, NN good lr'.format(s=size_large, nn=separation1),
                        edgelabel='Y',
                        color='blue',
                        qlabel='accept LR',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_nclustered_S_lr_NNlr.mask] = 1

    # compact not isolated, nnsmall, lr, NNnlr
    M_small_nisol_nclustered_S_lr_NNnlr = M_small_nisol_nclustered_S_lr.submask(lofarcat['NN_LR'] < lLR_thresh,
                        'NNnlr',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), good LR, NN bad lr'.format(s=size_large, nn=separation1),
                        edgelabel='N',
                        color='blue',
                        qlabel='accept LR',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_nclustered_S_lr_NNnlr.mask] = 1

    # compact not isolated, nnsmall, nlr
    M_small_nisol_nclustered_S_nlr = M_small_nisol_nclustered_S.submask(lofarcat['LR'] < lLR_thresh,
                        'nlr',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), bad LR'.format(s=size_large, nn=separation1),
                        edgelabel='N',
                        qlabel='NN LR > {l:.2f}?'.format(l=lLR_thresh),
                        masterlist=masterlist)

    # compact not isolated, nnsmall, nlr, NNlr
    M_small_nisol_nclustered_S_nlr_NNlr = M_small_nisol_nclustered_S_nlr.submask(lofarcat['NN_LR'] >= lLR_thresh,
                        'NNlr',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), bad LR, NN good lr'.format(s=size_large, nn=separation1),
                        edgelabel='Y',
                        color='red',
                        qlabel='accept no LR',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_nclustered_S_nlr_NNlr.mask] = 1

    # compact not isolated, nnsmall, nlr, NNnlr - there are possible doubles here!!
    M_small_nisol_nclustered_S_nlr_NNnlr = M_small_nisol_nclustered_S_nlr.submask(lofarcat['NN_LR'] < lLR_thresh,
                        'NNnlr',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), bad LR, NN bad lr'.format(s=size_large, nn=separation1),
                        edgelabel='N',
                        #color='orange',
                        qlabel='0.1 < flux ratio < 10?',
                        masterlist=masterlist)

    C1_simflux = (lofarcat['NN_Frat'] <= 10) & (lofarcat['NN_Frat'] > 0.1)
    M_small_nisol_nclustered_S_nlr_NNnlr_simflux = M_small_nisol_nclustered_S_nlr_NNnlr.submask(C1_simflux,
                        'simflux',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), bad LR, NN bad lr, sim flux'.format(s=size_large, nn=separation1),
                        edgelabel='Y',
                        #color='red',
                        qlabel='S1+S2 >= 50*(sep/100)**2 ?',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_nclustered_S_nlr_NNnlr_simflux.mask] = 1

    M_small_nisol_nclustered_S_nlr_NNnlr_diffflux = M_small_nisol_nclustered_S_nlr_NNnlr.submask(~C1_simflux,
                        'diffflux',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), bad LR, NN bad lr, diffflux'.format(s=size_large, nn=separation1),
                        edgelabel='N',
                        color='red',
                        qlabel='accept no LR',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_nclustered_S_nlr_NNnlr_diffflux.mask] = 1

    C2_dist = ((lofarcat['NN_Total_flux']+lofarcat['Total_flux']) >= 50*(lofarcat['NN_sep']/100.)**2.)
    M_small_nisol_nclustered_S_nlr_NNnlr_simflux_sep = M_small_nisol_nclustered_S_nlr_NNnlr_simflux.submask(C2_dist,
                        'dist',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), bad LR, NN bad lr, sim flux'.format(s=size_large, nn=separation1),
                        edgelabel='Y',
                        qlabel='visual confirmation?',
                        color='orange',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_nclustered_S_nlr_NNnlr_simflux_sep.mask] = 5
    
    if step == 4:
        M_small_nisol_nclustered_S_nlr_NNnlr_simflux_sep_lgz = M_small_nisol_nclustered_S_nlr_NNnlr_simflux_sep.submask(lofarcat['double_flag'] == 2,
                            'lgz',
                            edgelabel='complex',
                            color='green',
                            qlabel='lgz',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_small_nisol_nclustered_S_nlr_NNnlr_simflux_sep_lgz.mask] = 3
        
        M_small_nisol_nclustered_S_nlr_NNnlr_simflux_sep_lgz3 = M_small_nisol_nclustered_S_nlr_NNnlr_simflux_sep.submask(lofarcat['double_flag2'] == 2,
                            'lgz3',
                            edgelabel='complex',
                            color='green',
                            qlabel='lgz3',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_small_nisol_nclustered_S_nlr_NNnlr_simflux_sep_lgz3.mask] = 330
        
        
        M_small_nisol_nclustered_S_nlr_NNnlr_simflux_sep_noid = M_small_nisol_nclustered_S_nlr_NNnlr_simflux_sep.submask((lofarcat['double_flag'] == 1) | (lofarcat['double_flag2'] == 1),
                            'noid',
                            edgelabel='no match',
                            color='red',
                            qlabel='accept no LR',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_small_nisol_nclustered_S_nlr_NNnlr_simflux_sep_noid.mask] = 1
        
        
        M_small_nisol_nclustered_S_nlr_NNnlr_simflux_sep_art = M_small_nisol_nclustered_S_nlr_NNnlr_simflux_sep.submask((lofarcat['double_flag'] == 3) | (lofarcat['double_flag2'] == 3),
                            'artefact',
                            edgelabel='artefact',
                            color='grey',
                            qlabel='artefact',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_small_nisol_nclustered_S_nlr_NNnlr_simflux_sep_art.mask] = -1
        
        
        
        M_small_nisol_nclustered_S_nlr_NNnlr_simflux_sep_prob = M_small_nisol_nclustered_S_nlr_NNnlr_simflux_sep.submask(((lofarcat['double_flag'] == 4) | (lofarcat['double_flag'] == 0)) & (lofarcat['double_flag2'] == 0),
                            'prob',
                            edgelabel='problem - noid',
                            color='red',
                            qlabel='problem',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_small_nisol_nclustered_S_nlr_NNnlr_simflux_sep_prob.mask] = 1 ## these should be in lgz assoc! but let's double check...

    M_small_nisol_nclustered_S_nlr_NNnlr_simflux_nsep = M_small_nisol_nclustered_S_nlr_NNnlr_simflux.submask(~C2_dist,
                        'ndist',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), bad LR, NN bad lr, sim flux'.format(s=size_large, nn=separation1),
                        edgelabel='N',
                        color='red',
                        qlabel='accept no LR?',
                        masterlist=masterlist)    
    lofarcat['ID_flag'][M_small_nisol_nclustered_S_nlr_NNnlr_simflux_nsep.mask] = 1
    
    
    # other masks


    M_encloses = Mask(lofarcat['Encloses'],
                    'Encloses',
                    'encloses',
                    qlabel='encloses?')
    M_enclosed = Mask(lofarcat['Enclosed'],
                    'Enclosed',
                    'enclosed',
                    qlabel='enclosed?')
    M_intersects = Mask(lofarcat['Intersects'],
                    'Intersects',
                    'intersects',
                    qlabel='intersects?')
    M_encloses.make_sample(lofarcat)
    M_enclosed.make_sample(lofarcat)
    M_intersects.make_sample(lofarcat)
    

    #maskDC0 = lofarcat['DC_Maj'] == 0
    maskDC0 = lofarcat['Maj'] == 0

    M_S = Mask(lofarcat['S_Code'] == 'S', 'single')
    M_M = Mask(lofarcat['S_Code'] == 'M', 'multiple')
    M_C = Mask(lofarcat['S_Code'] == 'C', 'complex')

    M_Ngaus = []
    for i in range(1,6):
        M_Ngaus.append(Mask(lofarcat['Ng'] == i, 'Ng='+str(i)))
    M_Ngaus.append(Mask(lofarcat['Ng'] >= i+1, 'Ng>='+str(i+1)))

    M_huge = Mask(lofarcat['Maj'] > size_huge, 'huge')

    M_small = Mask(lofarcat['Maj'] <= size_large, 'small')

    M_isol = Mask(lofarcat['NN_sep'] > separation1, 'isol')

    M_cluster = Mask(lofarcat['NN5_sep'] < separation1, 'clustered',
                    'Clustered (5 sources within sep1)')

    M_bright = Mask(lofarcat['Total_flux'] > fluxcut, 'bright')
    M_nlr = Mask(lofarcat['LR'] >= lLR_thresh, 'lr')
    M_lr = Mask(lofarcat['LR'] < lLR_thresh,'nlr')

        
    M_huge = Mask(lofarcat['Maj'] > 100., 'huge')

    M_LGZ2 = Mask(lofarcat['LGZ_flag'] == 2., 'LGZv2')
                  
    #M_LGZz2 = Mask(lofarcat['LGZ_flag'] == 20., 'LGZv2_zoom')
                  

    # make a test sample for each final mask
    makesample = 1
    if makesample:
        for t in masterlist:
            if not t.has_children :
                print(t.name)
                t.make_sample(lofarcat,Nsample=None)

    # test that the final masks are indeed mutually disjoint and cover all sources
    endlist = []
    for t in masterlist:
        if not t.has_children:
            endlist.append(t)
    if not Masks_disjoint_complete(endlist):
        print('WARNING: children aren\'t disjoint and complete')


    # store the flags if we run in step 1 or stpe 2
    fcflg = 'FC_flag{st:d}'.format(st=step)
    if fcflg not in lofarcat.colnames:
        lofarcat.add_column(Column(-1*np.ones(len(lofarcat),dtype=int), fcflg))
    else:
        lofarcat[fcflg] = -1*np.ones(len(lofarcat),dtype=int)
        
        
    i = 0
    for t in masterlist:
        if not t.has_children:
            print(i, np.sum(lofarcat[fcflg][t.mask]==7))
            lofarcat[fcflg][t.mask] = i
            print(i, np.sum(t.mask))
            i += 1
            
    print(fcflg,'count')
    t,i = np.unique(lofarcat[fcflg], return_counts=True)
    for tt,ii in zip(t,i): print(tt,ii)
    

    ## write output file
    lofarcat.write(lofarcat_file_srt, overwrite=True)

    #sys.exit()

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

        A.add_node('start', label='ALL\n{n:n}'.format(n=M_all.N), shape='parallelogram') 
        #A.add_node('m_all', label='Large?'.format(n=M_all.N), shape='diamond') 

        A.add_edge('start', 'all', label='', penwidth=M_all.f*PW)
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
        A.draw('flow_step{st:d}_s{s:.0f}_nn{nn:.0f}.png'.format(st=step,s=size_large,nn=separation1)) # write to file
        A.draw('flow_step{st:d}_s{s:.0f}_nn{nn:.0f}.pdf'.format(st=step,s=size_large,nn=separation1)) # write to file
        A.write('flow_step{st:d}_s{s:.0f}_nn{nn:.0f}.dot'.format(st=step,s=size_large,nn=separation1)) # write to file




    print('ID_flag count')
    t,i = np.unique(lofarcat['ID_flag'], return_counts=True)
    for tt,ii in zip(t,i): print(tt,ii)
    
    print('ID_flag count - encloses')
    t,i = np.unique(lofarcat['ID_flag'][lofarcat['Encloses']], return_counts=True)
    for tt,ii in zip(t,i): print(tt,ii)
    print('ID_flag count - enclosed')
    t,i = np.unique(lofarcat['ID_flag'][lofarcat['Enclosed']], return_counts=True)
    for tt,ii in zip(t,i): print(tt,ii)
    print('ID_flag count - intersects')
    t,i = np.unique(lofarcat['ID_flag'][lofarcat['Intersects']], return_counts=True)
    for tt,ii in zip(t,i): print(tt,ii)
    
    
    print(fcflg,'count')
    t,i = np.unique(lofarcat[fcflg], return_counts=True)
    for tt,ii in zip(t,i): print(tt,ii)
    
    print('MC1_flag1 count')
    t,i = np.unique(lofarcat['MC1_flag1'], return_counts=True)
    for tt,ii in zip(t,i): print(tt,ii)
    print('MC2_flag1 count')
    t,i = np.unique(lofarcat['MC2_flag1'], return_counts=True)
    for tt,ii in zip(t,i): print(tt,ii)


    ## TESTING ##
    ## check gaus ML



if 0:
    fluxcuts = np.logspace(-4, 0, 1000)
    nS_fluxcuts = np.nan*np.zeros(len(fluxcuts))
    for fi,fluxcut in enumerate(fluxcuts):
        m = lofarcat['Total_flux']/1e3 > fluxcut
        nS_fluxcuts[fi] = 1.*np.sum(lofarcat['S_Code'][m] == 'S') /np.sum(m)
        #nS_fluxcuts[fi] = 1.*np.sum(m)
    f,ax = pp.paper_single_ax()
    ax.plot(fluxcuts, nS_fluxcuts)
    ax.set_ylabel('f(Single) ($S>S_{cut}$)')
    ax.set_xlabel('$\log S_{cut}$ [Jy]')
    plt.savefig('fraction_single_vs_S')


    sizecuts = np.linspace(15, 60, 10)
    fluxcuts = np.logspace(-3, 1, 1000)
    f,ax = pp.paper_single_ax()
    for si,sizecut in enumerate(sizecuts):
        ms = lofarcat['Maj'] > sizecut
        nS_fluxcuts = np.nan*np.zeros(len(fluxcuts))
        for fi,fluxcut in enumerate(fluxcuts):
            m = ms & (lofarcat['Total_flux']/1e3 > fluxcut)
            nS_fluxcuts[fi] = 1.*np.sum(m) /np.sum(ms)
        ax.plot(fluxcuts, nS_fluxcuts)
    #ax.set_ylabel('$f(Maj>Maj_{cut})$ ($S>S_{cut}$)')
    #ax.set_xlabel('$\log S_{cut}$ [Jy]')
    plt.savefig('fraction_large_vs_S')



    sizecuts = np.arange(10, 35, 1)
    NNcuts = np.arange(20, 125, 5)
    IM = np.zeros((len(sizecuts), len(NNcuts)))
    fluxcuts = np.logspace(-3, 1, 1000)
    f,ax = pp.paper_single_ax()
    for si,sizecut in enumerate(sizecuts):
        for ni,NNcut in enumerate(NNcuts):
            m = (lofarcat['Maj'] <= sizecut) & (lofarcat['NN_sep'] >= NNcut)
            IM[si,ni] = np.sum(m)
    IM = IM/Ncat
    c = ax.imshow(IM.T, origin='lower', extent=(10,60, 20,120))
    cbar = plt.colorbar(c)
    cbar.set_label('fraction')
    ax.invert_xaxis()
    ax.set_xlabel(r'$<$ size [arcsec]')
    ax.set_ylabel(r'$>$ NN separation [arcsec]')
    plt.savefig('number_compact_isolated')

    f,axs = plt.subplots(1,2,sharex=False,sharey=True,figsize=(12,6))
    ax=axs[0]
    ax.plot(NNcuts,IM.T)
    ax.set_ylabel('fraction')
    ax.set_xlabel(r'$>$ NN separation [arcsec]')
    ax=axs[1]
    ax.plot(sizecuts,IM)
    ax.set_xlabel(r'$<$ size [arcsec]')


    nb=100
    if step == 2:
        # plot LR distribuion for different classes
        f,ax = pp.paper_single_ax()
        _ =ax.hist(np.log10(1.+lofarcat['LR']), bins=100, normed=True, log=False,histtype='step',color='k',linewidth=2,label='All')
        _ =ax.hist(np.log10(1.+lofarcat['LR'][M_small_isol_S.mask]), bins=100, normed=True, histtype='step', label=M_small_isol_S.name.replace('_','\_'))
        #_ =ax.hist(np.log10(1.+lofarcat['LR'][m_small_isol_nS]), bins=100, normed=True, histtype='step', label=l_small_isol_nS)
        _ =ax.hist(np.log10(1.+lofarcat['LR'][M_small_nisol.mask]), bins=100, normed=True, histtype='step', label=M_small_nisol.name.replace('_','\_'))
        _ =ax.hist(np.log10(1.+lofarcat['LR'][M_large.mask]), bins=100, normed=True, histtype='step', label=M_large.name.replace('_','\_'))
        _ =ax.hist(np.log10(1.+lofarcat['LR'][M_large_faint_complex.mask]), bins=100, normed=True, histtype='step', label=M_large_faint_complex.name.replace('_','\_'))
        #_ =ax.hist(np.log10(1.+lofarcat['LR'][M_large_faint_nhuge_n2masx.mask]), bins=100, normed=True, histtype='step', label=M_large_faint_nhuge_n2masx.name.replace('_','\_'))
        ax.legend()
        ax.set_ylim(0,2)
        ax.set_xlabel('$\log (1+LR)$')
        ax.set_ylabel('$N$')
        plt.savefig('lr_dist_classes')

    # plot LR distribuion for different classes
    f,ax = pp.paper_single_ax()
    counts, xedges, yedges, im =ax.hist2d(np.log10(1.+lofarcat['LR']), np.log10(lofarcat['Maj']), bins=100, density=True, vmin=0, vmax=2, label='')
    cbar = plt.colorbar(im, ax=ax)
    ax.legend()
    #ax.set_ylim(0,2)
    ax.set_xlabel('$\log (1+LR)$')
    ax.set_ylabel('$\log$ Maj [arcsec]')
    cbar.set_label('$N$')
    plt.savefig('lr_dist_size')

    f,ax = pp.paper_single_ax()
    counts, xedges, yedges, im =ax.hist2d(np.log10(1.+lofarcat['LR']), np.log10(lofarcat['Total_flux']), bins=100, density=True, vmin=0, vmax=2, label='')
    cbar = plt.colorbar(im, ax=ax)
    ax.legend()
    #ax.set_ylim(0,2)
    ax.set_xlabel('$\log (1+LR)$')
    ax.set_ylabel('$\log S$ [Jy]')
    cbar.set_label('$N$')
    plt.savefig('lr_dist_flux')



    #f,ax = pp.paper_single_ax()
    #f = plt.figure()
    f,axs = plt.subplots(1,2,sharex=True,sharey=True,figsize=(12,6))
    ax = axs[0]
    counts, xedges, yedges, im =ax.hist2d(np.log10(lofarcat['Maj'][M_S.mask]), np.log10(lofarcat['Total_flux'][m_S]), bins=100, label='')
    cbar = plt.colorbar(im, ax=ax)
    x1,x2 = ax.get_xlim()
    y1,y2 = ax.get_ylim()
    ax.vlines(np.log10(15.),y1,y2)
    ax.hlines(np.log10(10.),x1,x2)
    ax.legend()
    ax.set_title('S')
    #ax.set_ylim(0,2)
    ax.set_xlabel('$\log $ Maj [arcsec]')
    ax.set_ylabel('$\log S$ [mJy]')
    cbar.set_label('$N$')
    ax = axs[1]
    counts, xedges, yedges, im =ax.hist2d(np.log10(lofarcat['Maj'][~M_S.mask]), np.log10(lofarcat['Total_flux'][~m_S]), bins=100, label='')
    cbar = plt.colorbar(im, ax=ax)
    x1,x2 = ax.get_xlim()
    y1,y2 = ax.get_ylim()
    ax.vlines(np.log10(15.),y1,y2)
    ax.hlines(np.log10(10.),x1,x2)
    ax.legend()
    ax.set_title('!S')
    #ax.set_ylim(0,2)
    ax.set_xlabel('$\log $ Maj [arcsec]')
    #ax.set_ylabel('$\log S$ [mJy]')
    cbar.set_label('$N$')
    plt.savefig('lr_dist_size_flux')




    #f,ax = pp.paper_single_ax()
    #f = plt.figure()
    f,axs = plt.subplots(1,2,sharex=True,sharey=True,figsize=(12,6))
    ax = axs[0]
    counts, xedges, yedges, im =ax.hist2d((lofarcat['Maj'][M_lr.mask]), (lofarcat['NN_sep'][M_lr.mask]), bins=200, range=((0,50),(0,200)), label='')
    cbar = plt.colorbar(im, ax=ax)
    x1,x2 = ax.get_xlim()
    y1,y2 = ax.get_ylim()
    ax.vlines((15.),y1,y2)
    ax.hlines((10.),x1,x2)
    ax.legend()
    ax.set_title('good LR')
    #ax.set_ylim(0,2)
    ax.set_xlabel('Maj [arcsec]')
    ax.set_ylabel('NN separation [arcsec]')
    cbar.set_label('$N$')
    ax = axs[1]
    counts, xedges, yedges, im =ax.hist2d((lofarcat['Maj'][~M_lr.mask]), (lofarcat['NN_sep'][~M_lr.mask]), bins=200, range=((0,50),(0,200)), label='')
    cbar = plt.colorbar(im, ax=ax)
    x1,x2 = ax.get_xlim()
    y1,y2 = ax.get_ylim()
    ax.vlines((15.),y1,y2)
    ax.hlines((10.),x1,x2)
    ax.legend()
    ax.set_title('bad LR')
    #ax.set_ylim(0,2)
    ax.set_xlabel('Maj [arcsec]')
    #ax.set_ylabel('$\log S$ [mJy]')
    cbar.set_label('$N$')
    plt.savefig('lr_dist_size_nnsep')





    # # diagnostic plots 


    # plot size distribution
    f, ax = pp.paper_single_ax()
    ax.hist(lofarcat['Maj'][~maskDC0], range=(0,80), bins=100, histtype='step', label='All')
    ax.hist(lofarcat['Maj'][~maskDC0&M_S.mask], range=(0,80), bins=100, histtype='step', label='S')
    ax.hist(lofarcat['Maj'][~maskDC0&M_M.mask], range=(0,80), bins=100, histtype='step', label='M')
    ax.hist(lofarcat['Maj'][~maskDC0&M_C.mask], range=(0,80), bins=100, histtype='step', label='C')
    ax.set_xlabel('Major Axis [arcsec]')
    ax.set_ylabel('N')
    ax.legend()
    plt.savefig('size_dist_classes')


    # plot nearest neighbour distribution
    f,ax = pp.paper_single_ax()
    ax.hist(lofarcat['NN_sep'], bins=100, histtype='step', label='All')
    ax.hist(lofarcat['NN_sep'][M_S.mask], bins=100, histtype='step', label='S')
    ax.set_xlabel('Nearest source [arcsec]')
    ax.set_ylabel('N')
    ax.legend()
    plt.savefig('NNdist_dist')


    # 2D histogram : size-nearest neighbour distance
    # for 'S' sources
    f,ax = pp.paper_single_ax()
    X =  lofarcat['NN_sep'][~maskDC0&M_S.mask]
    Y = lofarcat['Maj'][~maskDC0&M_S.mask]
    H, xe, ye =  np.histogram2d( X, Y, bins=(100,100), normed=True)
    H2 = H.T
    xc = (xe[1:] +xe[:-1] )/2.
    yc = (ye[1:] +ye[:-1] )/2.
    c = ax.contour(xc, yc, H2, [0.5])
    xind = np.sum(X>xe[:,np.newaxis],axis=0)-1
    yind = np.sum(Y>ye[:,np.newaxis],axis=0)-1
    Hval = H2[yind,xind]
    c = ax.scatter(X, Y,c=Hval,s=10, edgecolor='none',zorder=1)
    x1,x2 = ax.get_xlim()
    y1,y2 = ax.get_ylim()
    ax.hlines(size_large,x1,x2,colors='k',linestyle='dashed')
    ax.vlines(separation1,y1,y2,colors='k',linestyle='dashed')
    ax.set_xlabel('NN separation [arcsec]')
    ax.set_ylabel('DCmaj [arcsec]')
    ax.contour(xc, yc, H2)
    plt.savefig('size_NNdist_dist_s')


    # and 'M' sources
    f,ax = pp.paper_single_ax()
    X =  lofarcat['NN_sep'][~maskDC0&M_M.mask]
    Y = lofarcat['Maj'][~maskDC0&M_M.mask]
    H, xe, ye =  np.histogram2d( X, Y, bins=(100,100), normed=True)
    H2 = H.T
    xc = (xe[1:] +xe[:-1] )/2.
    yc = (ye[1:] +ye[:-1] )/2.
    c = ax.contour(xc, yc, H2, [0.5])
    xind = np.sum(X>xe[:,np.newaxis],axis=0)-1
    yind = np.sum(Y>ye[:,np.newaxis],axis=0)-1
    Hval = H2[yind,xind]
    c = ax.scatter(X, Y,c=Hval,s=10, edgecolor='none',zorder=1)
    x1,x2 = ax.get_xlim()
    y1,y2 = ax.get_ylim()
    ax.hlines(size_large,x1,x2,colors='k',linestyle='dashed')
    ax.vlines(separation1,y1,y2,colors='k',linestyle='dashed')
    ax.set_xlabel('NN separation [arcsec]')
    ax.set_ylabel('DCmaj [arcsec]')
    ax.contour(xc, yc, H2)
    plt.savefig('size_NNdist_dist_m')



