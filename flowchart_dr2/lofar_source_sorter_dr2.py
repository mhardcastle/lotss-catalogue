'''
lofar_source_sorter.py
the main code for the decision tree
'''
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from shapely.geometry import Polygon, Point

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


def inHETDEX(ra,dec):
    inside = np.zeros(len(ra), dtype=bool)
    
    hvert = Table.read('/home/wwilliams/data2/projects/lofar_surveys/DR1/hetdex_vertices.fits')
    hra = hvert['HETDEX_RA'][0]
    hdec = hvert['HETDEX_DEC'][0]
    
    phetdex = Polygon(np.array((hra, hdec)).T)
    for i in range(len(ra)):
        inside[i] = phetdex.contains(Point(ra[i],dec[i]))
    
    return inside

if __name__=='__main__':

    if len(sys.argv) == 1:
        print("Usage is : python lofar_source_sorter_dr2.py field_code step_no ")
        print('E.g.: python lofar_source_sorter_dr2.py 0 1 ')
        sys.exit(1)

    h = str(sys.argv[1])
    if 'h' not in h:
        h+='h'
    if h not in  ['0h','13h']:
        print('unknown field code (should be 0h or 13h)',h)
        sys.exit(1)
        
    step = int(sys.argv[2])
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
    lLR_thresh = 0.309            # LR threshold
    fluxcut = 8.               # in mJy
    #fluxcut2 = 4.0               # in mJy
    #fluxcut2 = 3.0               # in mJy
    #fluxcut2 = 5.0               # in mJy
    fluxcut2 = float(sys.argv[3])
    
    ##LR threshold
    ##s0h -> 0.394
    ##s13h -> 0.328
    ##n13h -> 0.309 

    ### Required INPUTS
    # lofar source catalogue, gaussian catalogue and ML catalogues for each

    path = '/data2/wwilliams/projects/lofar_surveys/LoTSS-DR2-Feb2020/'
    
    version ='v100'
    
    if step == 1:
        # lr infor is in the catalogue - this makes things easier
        lofarcat_file = path+'LoTSS_DR2_{version}.srl_{h}.lr-full.presort.hdf5'.format(h=h,version=version)
        #lofargcat_file = path+'lr/LoTSS_DR2_{version}.gaus_{h}.lr-full.fits'.format(h=h,version=version)
    elif step >= 2:
        lofarcat_file = path+'LoTSS_DR2_{version}.srl_{h}.lr-full.sorted_step1_flux{ff:.0f}.hdf5'.format(h=h,version=version,ff=fluxcut2)
        #lofargcat_file = path+'lr/LoTSS_DR2_{version}.gaus_{h}.lr-full.fits'.format(h=h,version=version)
        


    # output file - with updated columns
    lofarcat_file_srt = path+'LoTSS_DR2_{version}.srl_{h}.lr-full.sorted_step{st}_flux{ff:.0f}.hdf5'.format(h=h,version=version,st=step,ff=fluxcut2)



    ## Gaus catalogue
    #lofargcat = Table.read(lofargcat_file)
    # only relevant gaussians are in M or C sources
    #lofargcat = lofargcat[lofargcat['S_Code'] != 'S']

    # Source catalogue
    lofarcat = Table.read(lofarcat_file)
    
    if 'RA_1' in lofarcat.colnames:
        lofarcat.rename_column('RA_1','RA')
        lofarcat.rename_column('DEC_1','DEC')
        lofarcat.rename_column('LR_1a','LR')
    
    Nlofarcat = len(lofarcat)
    
    
    ## write output file
    #lofarcat.write(lofarcat_file_srt, overwrite=True)
    #sys.exit()

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
        if 'ML_flag' not in lofarcat.colnames:
            # add a temporary ML_flag - set all to 0 for LGZ
            lofarcat.add_column(Column(data=np.zeros(len(lofarcat),dtype=bool), name='ML_flag'))
            print('added ML flag to send everything to LGZ temporarily')
    elif step == 2:
        print('msources have been handled, now to give them the right ID flags')
        if 'msource_flag1' not in lofarcat.colnames:
            raise  RuntimeError('need the msource_flag1 -- run get_msource_flag.py')
        
        #fixlist = ['ILTJ111611.15+493234.8', 'ILTJ121557.43+512418.5', 'ILTJ121956.33+473756.2', 'ILTJ124047.70+500956.4', 'ILTJ131951.11+531748.0', 'ILTJ135441.11+541646.4']
        #for s in fixlist:
            #lofarcat['msource1_flag'][lofarcat['Source_Name']==s] = 1
        #fixlist = ['ILTJ140430.41+560657.9'] 
        #for s in fixlist:
            #lofarcat['msource1_flag'][lofarcat['Source_Name']==s] = 2
    elif step == 3:
        print('need ML outputs')  ### actually, these are now there from the start
        print('machine learning has been run, now to give them the right ID flags')
        if 'ML_flag' not in lofarcat.colnames:
            raise  RuntimeError('need the ML_flag -- run get_ml_flags.py')  ## 1 for LR, 0 for LGZ
    else:
        print('step {s} not defined, quitting'.format(s=step))
        sys.exit(1)

    if step == 1:
        lofarcat.add_column(Column(data=np.zeros(len(lofarcat),dtype=bool), name='WEAVE_priority1'))
        lofarcat.add_column(Column(data=np.zeros(len(lofarcat),dtype=bool), name='WEAVE_priority2'))
        lofarcat.add_column(Column(data=np.zeros(len(lofarcat),dtype=bool), name='WEAVE_priority3'))
        ind8hr = (lofarcat['RA']>= 110) & (lofarcat['RA']<= 135) & (lofarcat['DEC']>= 26) & (lofarcat['DEC']<= 41)
        ind13hr45 = (lofarcat['RA']>= 130) & (lofarcat['RA']<= 280) & (lofarcat['DEC']>= 40) & (lofarcat['DEC']<= 45)
        ind13hr60 = (lofarcat['RA']>= 120) & (lofarcat['RA']<= 240) & (lofarcat['DEC']>= 60) & (lofarcat['DEC']<= 65)
        
        #dr1_mos_ids = ['P10Hetdex', 'P11Hetdex12', 'P12Hetdex11', 'P14Hetdex04', 'P15Hetdex13', 'P164+55', 'P169+55', 'P16Hetdex13', 'P173+55', 'P178+55', 'P182+55', 'P187+55', 'P18Hetdex03', 'P191+55', 'P196+55', 'P19Hetdex17', 'P1Hetdex15', 'P200+55', 'P205+55', 'P206+50', 'P206+52', 'P209+55', 'P21', 'P210+47', 'P211+50', 'P213+47', 'P214+55', 'P217+47', 'P218+55', 'P219+50', 'P219+52', 'P221+47', 'P223+50', 'P223+52', 'P223+55', 'P225+47', 'P227+50', 'P227+53', 'P22Hetdex04', 'P23Hetdex20', 'P25Hetdex09', 'P26Hetdex03', 'P27Hetdex09', 'P29Hetdex19', 'P30Hetdex06', 'P33Hetdex08', 'P34Hetdex06', 'P35Hetdex10', 'P37Hetdex15', 'P38Hetdex07', 'P39Hetdex19', 'P3Hetdex16', 'P41Hetdex', 'P42Hetdex07', 'P4Hetdex16', 'P6', 'P7Hetdex11', 'P8Hetdex']
        #indhetdex = np.zeros(len(lofarcat),dtype=bool)
        #for dr1_m in dr1_mos_ids:
            #indm = lofarcat['Mosaic_ID'] == dr1_m
            #indhetdex[indm] = True
            
            
        
        lofarcat['WEAVE_priority1'][ind13hr60] = True
        
        print('{} sources out of {} with WEAVE_priority1 ({:.1f}%)'.format(np.sum(lofarcat['WEAVE_priority1']), Nlofarcat, 100.*np.sum(lofarcat['WEAVE_priority1'])/Nlofarcat))
        
        
        lofarcat['WEAVE_priority3'][ind8hr] = True  
        lofarcat['WEAVE_priority2'][ind13hr45] = True
        
        ## remove dr1 area:   - temp take out - slow and small area
        #indhetdex = inHETDEX(lofarcat['RA'],lofarcat['DEC'])
        #lofarcat['WEAVE_priority1'][indhetdex] = False
        print('{} sources out of {} with WEAVE_priority1 not hetdex ({:.1f}%)'.format(np.sum(lofarcat['WEAVE_priority1']), Nlofarcat, 100.*np.sum(lofarcat['WEAVE_priority1'])/Nlofarcat))


    # this is easy to run...
    # get the large 2masx sources (must run match_2masx for these)
    if '2MASX_match_large' not in lofarcat.colnames:
        raise  RuntimeError('need the 2masx information -- run match_2masx.py')
    big2masx = lofarcat['2MASX_match_large']

    # step 1 happens before we have visual classifications - if any ever for dr2
    if step in [1,2,3]:
        # we start with no flagging - visual, other
        cleanflag = np.ones(Nlofarcat,dtype=bool)
        flag0 = np.zeros(len(lofarcat),dtype=bool)
        
        clustered_flag = flag0
        Lclustered_flag = flag0
        huge_faint_flag = flag0
        nhuge_faint_flag = flag0
        nhuge_2masx_flag = flag0
        double_flag = flag0
        double_flag2 = flag0
        mnisol_flag2 = flag0
        mnisol_flag1 = flag0
        #big2masx = flag0
        edge_flag = flag0
        
        if 'artefact_flag' in lofarcat.colnames:
            artefact_flag = lofarcat['artefact_flag']
        else:
            artefact_flag = flag0
        
        # combine the artefact flags
        # artefacts have been identified through various routes of visual checking
        Artefact_flag = (artefact_flag == 1)
        
        if '2MASX' not in lofarcat.colnames:
            lofarcat.add_column(Column(data=flag0, name='2MASX'))
        if 'double_flag' not in lofarcat.colnames:
            lofarcat.add_column(Column(data=flag0, name='double_flag'))
        
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

    # 1
    M_all_full = Mask(lofarcat['RA'] > -1,
                    'all',
                    qlabel='weave priority?',
                    masterlist=masterlist,
                    maskflag=lofarcat['ML_flag'], maskflagstr='LR:LGZ')


    
    # non-weave
    # 10
    M_all = M_all_full.submask(lofarcat['WEAVE_priority1']==True,
                        'weave',
                        'Weave priority',
                        edgelabel='Y',
                        color='gray',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_all.mask] = -1


    #non -weave
    # 11
    M_all_noweave = M_all_full.submask(lofarcat['WEAVE_priority1']==False,
                        'no_weave',
                        'No Weave priority',
                        edgelabel='N',
                        color='gray',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_all_noweave.mask] = -1


    # sources 
    # 100
    M_all_clean = M_all.submask(~artefact_flag,
                        'src',
                        'Clean',
                        edgelabel='N',
                        qlabel='Huge 2MASX?\n(r(2MASX)>60")',
                        masterlist=masterlist)


    # artefacts
    # 101
    M_all_artefact = M_all.submask(artefact_flag,
                        'artefact',
                        'Artefact\n(visually confirmed)',
                        edgelabel='Y',
                        color='gray',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_all_artefact.mask] = -1


    # sources - clean
    # 1000
    M_all_clean2 = M_all_clean.submask(~big2masx,
                        'clean',
                        'Clean2',
                        edgelabel='N',
                        qlabel='Large?\n(s>{s:.0f}")'.format(s=size_large),
                        masterlist=masterlist)



    # big optical gal 
    # 1001
    M_all_biggal = M_all_clean.submask((big2masx),
                        'big2MASX',
                        'Huge 2MASX source\n(64 r(2MASX)>60" galaxies)',
                        edgelabel='Y',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_all_biggal.mask] = 2
    

    # large 
    # 10000
    M_large = M_all_clean2.submask(lofarcat['Maj'] > size_large,
                        'large',
                        'large (s>{s:.0f}")'.format(s=size_large),
                        qlabel='Bright?\n(S>{f:.0f} mJy)'.format(f=fluxcut, s=size_large),
                        masterlist=masterlist)

    # large bright
    # 100000
    M_large_bright = M_large.submask(lofarcat['Total_flux'] >= fluxcut,
                        'bright',
                        'large (s>{s:.0f}") & bright (S>={f:.0f} mJy)'.format(f=fluxcut, s=size_large),
                        qlabel='LGZ',
                        color='green',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_large_bright.mask] = 31
    lofarcat['LGZ_flag'][M_large_bright.mask] = 1

    # large faint
    # 100001
    M_large_faint = M_large.submask(lofarcat['Total_flux'] < fluxcut,
                        'faint',
                        'large (s>{s:.0f}") & faint (S<{f:.0f} mJy)'.format(f=fluxcut, s=size_large),
                        edgelabel='N',
                        qlabel='Brightish?\n(S>{f2:.0f} mJy)'.format(f2=fluxcut2, s=size_large),
                        masterlist=masterlist)
    
    lofarcat['ID_flag'][M_large_faint.mask] = 5  #5 is tbd

    # large faint
    # 1000010
    M_large_faint_faint = M_large_faint.submask(lofarcat['Total_flux'] < fluxcut2,
                        'faintfaint',
                        'large (s>{s:.0f}") & faint faint (S<{f2:.0f} mJy)'.format(f2=fluxcut2, s=size_large),
                        edgelabel='N',
                        qlabel='TBD',
                        color='orange',
                        masterlist=masterlist)
    
    lofarcat['ID_flag'][M_large_faint_faint.mask] = 5  #5 is tbd


    # large faint
    # 1000011
    M_large_mid_faint = M_large_faint.submask(lofarcat['Total_flux'] >= fluxcut2,
                        'midfaint',
                        'large (s>{s:.0f}") & mid faint ({f2:.0f}<=S<{f:.0f} mJy)'.format(f=fluxcut, f2=fluxcut2, s=size_large),
                        edgelabel='Y',
                        qlabel='ML output LGZ?',
                        #color='cyan',
                        masterlist=masterlist)
    
    lofarcat['ID_flag'][M_large_mid_faint.mask] = 5  #5 is tbd



    # large faint
    # 10000110
    M_large_mid_faint_mllgz = M_large_mid_faint.submask((lofarcat['ML_flag'] ==False),
                        'mllgz',
                        'ml lgz',
                        edgelabel='Y',
                        qlabel='LGZ',
                        color='green',
                        masterlist=masterlist)
    
    lofarcat['ID_flag'][M_large_mid_faint_mllgz.mask] = 3  #5 is tbd
    
    
    # large faint
    # 10000111
    M_large_mid_faint_nmllgz = M_large_mid_faint.submask((lofarcat['ML_flag'] ==True),
                        'nmllgz',
                        'no ml lgz',
                        edgelabel='N',
                        qlabel='Visual sorting',
                        color='cyan',
                        masterlist=masterlist)
    
    lofarcat['ID_flag'][M_large_mid_faint_nmllgz.mask] = 5  #5 is tbd




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
                            color='blue',
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
    # 10001
    M_small = M_all_clean2.submask(lofarcat['Maj'] <= size_large,
                        'small',
                        'compact (s<{s:.0f}")'.format(s=size_large),
                        edgelabel='N',
                        qlabel='Isolated?\n(NN>{nn:.0f}")'.format(nn=separation1),
                        masterlist=masterlist)
    # compact isolated
    # 100010
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
                        color='blue',
                        qlabel='Accept no LR',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_isol_S_nlr.mask] = 1


    # compact isolated nS
    M_small_isol_nS = M_small_isol.submask(lofarcat['S_Code'] != 'S',
                        'nS',
                        'compact isolated (s<{s:.0f}", NN>{nn:.0f}") !S'.format(s=size_large, nn=separation1),
                        edgelabel='N',
                        qlabel='Bright?\n(S>{f:.0f} mJy)'.format(f=fluxcut, s=size_large),
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_isol_nS.mask] = 5  # - these should all be overwritten below, but leave it here to doublecheck (there should be no 5's in the end)



    if step==1:
        # compact isolated nS
        M_small_isol_nS_bright = M_small_isol_nS.submask(lofarcat['Total_flux'] >= fluxcut,
                            'nS_bright',
                            'compact isolated bright (s<{s:.0f}", NN>{nn:.0f}", S> {f:.0f}) !S'.format(s=size_large, nn=separation1,f=fluxcut),
                            edgelabel='Y',
                            color='gray',
                            qlabel='msource',
                            masterlist=masterlist)
        M_small_isol_nS_faint = M_small_isol_nS.submask(lofarcat['Total_flux'] < fluxcut,
                            'nS_faint',
                            'compact isolated faint (s<{s:.0f}", NN>{nn:.0f}", S < {f:.0f}) !S'.format(s=size_large, nn=separation1,f=fluxcut),
                            edgelabel='N',
                            color='orange',
                            qlabel='msource TBD',
                            masterlist=masterlist)
    


    if step == 2:
        
        # compact isolated nS
        M_small_isol_nS_bright = M_small_isol_nS.submask(lofarcat['Total_flux'] >= fluxcut,
                            'nS_bright',
                            'compact isolated bright (s<{s:.0f}", NN>{nn:.0f}", S> {f:.0f}) !S'.format(s=size_large, nn=separation1,f=fluxcut),
                            edgelabel='Y',
                            qlabel='msource outcome?',
                            masterlist=masterlist)
        
        
        M_small_isol_nS_faint = M_small_isol_nS.submask(lofarcat['Total_flux'] < fluxcut,
                            'nS_faint',
                            'compact isolated faint (s<{s:.0f}", NN>{nn:.0f}", S < {f:.0f}) !S'.format(s=size_large, nn=separation1,f=fluxcut),
                            edgelabel='N',
                            color='orange',
                            qlabel='msource TBD',
                            masterlist=masterlist)
        
        
        ### msource1_flag
        #0: no match
        #1: accept ML of the source
        #2: accept ML of the gaussian with highest ML
        #3: deblend and accept both gaussians
        #4: deblend workflow
        #5: LOFAR galaxy zoo 
        #6: visual check 
        
        M_small_isol_nS_bright_LR = M_small_isol_nS_bright.submask((lofarcat['msource_flag1'] == 0) |
                                       (lofarcat['msource_flag1'] == 1) |
                                       (lofarcat['msource_flag1'] == 2),
                            'id_by_lr',
                            'id or not by LR',
                            edgelabel='1',
                            color='blue',
                            qlabel='accept LR id/no id',
                            masterlist=masterlist)
        
        M_small_isol_nS_bright_deblend = M_small_isol_nS_bright.submask((lofarcat['msource_flag1'] == 3) |
                                       (lofarcat['msource_flag1'] == 4),
                            'deblend',
                            'deblend',
                            edgelabel='2',
                            color='yellow',
                            qlabel='deblend',
                            masterlist=masterlist)
        
        M_small_isol_nS_bright_lgz = M_small_isol_nS_bright.submask((lofarcat['msource_flag1'] == 5),
                            'lgz',
                            'lgz',
                            edgelabel='3',
                            color='green',
                            qlabel='lgz',
                            masterlist=masterlist)
        
        M_small_isol_nS_bright_prefilt = M_small_isol_nS_bright.submask((lofarcat['msource_flag1'] == 6),
                            'prefilt',
                            'prefilt',
                            edgelabel='4',
                            color='cyan',
                            qlabel='prefilt',
                            masterlist=masterlist)
        
        
        ## TODO fix flags
        # the ID_flags here depend on the output of Philips flow chart:
        lofarcat['ID_flag'][M_small_isol_nS_bright_LR.mask] = 12
        lofarcat['ID_flag'][M_small_isol_nS_bright_deblend.mask] = 41
        lofarcat['ID_flag'][M_small_isol_nS_bright_lgz.mask] = 32
        lofarcat['ID_flag'][M_small_isol_nS_bright_prefilt.mask] = 32
        lofarcat['LGZ_flag'][M_small_isol_nS_bright.mask] = 20


    # compact not isolated
    # 100011
    M_small_nisol = M_small.submask(lofarcat['NN_sep'] <= separation1,
                        'nisol',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}")'.format(s=size_large, nn=separation1),
                        edgelabel='N',
                        qlabel='S?',
                        masterlist=masterlist)
    
    
    # compact not isolated, S
    # 1000110
    M_small_nisol_S = M_small_nisol.submask(lofarcat['S_Code'] == 'S',
                        'S',
                        edgelabel='Y',
                        #qlabel='Brightish?\n(S>{f2:.0f} mJy)'.format(f2=fluxcut2),
                        qlabel='Clustered?\n(NN5<{nn:.0f}")'.format(s=size_large, nn=separation1),
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_S.mask] = 5  # - these should all be overwritten below, but leave it here to doublecheck (there should be no 5's in the end)

    # compact not isolated, not S
    # 1000111
    M_small_nisol_nS = M_small_nisol.submask(lofarcat['S_Code'] != 'S',
                        'nS',
                        edgelabel='N',
                        qlabel='Brightish?\n(S>{f2:.0f} mJy)'.format(f2=fluxcut2),
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_nS.mask] = 5  # - these should all be overwritten below, but leave it here to doublecheck (there should be no 5's in the end)

    
    # compact not isolated, S, clustered
    M_small_nisol_S_clustered = M_small_nisol_S.submask((lofarcat['NN5_sep'] <= separation1),
                        'clustered',
                        'compact not isolated S clustered (s<{s:.0f}", NN5>{nn:.0f}")'.format(s=size_large, nn=separation1),
                        edgelabel='Y',
                        qlabel='Brightish?\n(S>{f2:.0f} mJy)'.format(f2=fluxcut2),
                        masterlist=masterlist)
    
    
    # compact not isolated, S, not clustered
    M_small_nisol_S_nclustered = M_small_nisol_S.submask((lofarcat['NN5_sep'] > separation1),
                        'nclustered',
                        'compact not isolated S not clustered (s<{s:.0f}", NN5>{nn:.0f}")'.format(s=size_large, nn=separation1),
                        edgelabel='N',
                        qlabel='LR > {l:.2f}?'.format(l=lLR_thresh),
                        masterlist=masterlist)
    
    

    # brightish
    M_small_nisol_S_clustered_faint = M_small_nisol_S_clustered.submask((lofarcat['Total_flux'] <= fluxcut2),
                        'tbc_faint',
                        edgelabel='N',
                        qlabel='TBD',
                        color='orange',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_S_clustered_faint.mask] = 5
    
    
    # faint
    M_small_nisol_S_clustered_brightish = M_small_nisol_S_clustered.submask((lofarcat['Total_flux'] >  fluxcut2),
                        'tbc_brightish',
                        edgelabel='Y',
                        qlabel='ML output LGZ?',
                        #qlabel='TBC',
                        #color='cyan',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_S_clustered_brightish.mask] = 5
    
    
    # faint -ML LGZ
    M_small_nisol_S_clustered_brightish_mllgz = M_small_nisol_S_clustered_brightish.submask((lofarcat['ML_flag'] ==False),
                        'ml_lgz',
                        edgelabel='Y',
                        qlabel='LGZ',
                        color='green',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_S_clustered_brightish_mllgz.mask] = 3
    
    
    # faint - ML not LGZ
    M_small_nisol_S_clustered_brightish_nmllgz = M_small_nisol_S_clustered_brightish.submask((lofarcat['ML_flag'] ==True),
                        'not_ml_lgz',
                        edgelabel='N',
                        qlabel='Visual sorting',
                        color='cyan',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_S_clustered_brightish_nmllgz.mask] = 5
        
        
        
    if step == 4:
        
        ## expand the M_small_nisol_S_clustered_faint based on visual flags
        pass

    # compact not isolated, nnlarge
    M_small_nisol_nS_bright = M_small_nisol_nS.submask(lofarcat['Total_flux'] >= fluxcut2,
                        'nS_brightish',
                        edgelabel='Y',
                        qlabel='msource',
                        color='gray',
                        masterlist=masterlist)
                        
    # compact not isolated, nnlarge
    M_small_nisol_nS_faint = M_small_nisol_nS.submask(lofarcat['Total_flux'] < fluxcut2,
                        'nS_faint',
                        edgelabel='N',
                        qlabel='msource TBD',
                        color='orange',
                        masterlist=masterlist)
    
    if step >= 2:
        
        # compact not isolated, nnlarge
        M_small_nisol_nS_bright_LR = M_small_nisol_nS_bright.submask((lofarcat['msource_flag1'] == 0) |
                                                   (lofarcat['msource_flag1'] == 1) |
                                                   (lofarcat['msource_flag1'] == 2), 
                            'lr',
                            edgelabel='lr',
                            qlabel='accept LR id/no id',
                            color='blue',
                            masterlist=masterlist)
        
        M_small_nisol_nS_bright_deblend = M_small_nisol_nS_bright.submask((lofarcat['msource_flag1'] == 3) |
                                                   (lofarcat['msource_flag1'] == 4), 
                            'deblend',
                            edgelabel='deblend',
                            qlabel='deblend',
                            color='yellow',
                            masterlist=masterlist)
        
        M_small_nisol_nS_bright_lgz = M_small_nisol_nS_bright.submask((lofarcat['msource_flag1'] == 5), 
                            'lgz',
                            edgelabel='lgz',
                            qlabel='lgz',
                            color='green',
                            masterlist=masterlist)
        
        M_small_nisol_nS_bright_prefilt = M_small_nisol_nS_bright.submask((lofarcat['msource_flag1'] == 6), 
                            'prefilt',
                            edgelabel='prefilt',
                            qlabel='prefilt',
                            color='cyan',
                            masterlist=masterlist)
        
        # TODO fix ID_flags
        # the ID_flags here depend on the output of Philips flow chart:
        lofarcat['ID_flag'][M_small_nisol_nS_bright_prefilt.mask ] = 12
        lofarcat['LGZ_flag'][M_small_nisol_nS_bright_prefilt.mask] = 20

    
    
    

    # compact not isolated, nnsmall, nlr
    M_small_nisol_S_nclustered_nlr = M_small_nisol_S_nclustered.submask(lofarcat['LR'] < lLR_thresh,
                        'nlr',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), bad LR'.format(s=size_large, nn=separation1),
                        edgelabel='N',
                        qlabel='NN LR > {l:.2f}?'.format(l=lLR_thresh),
                        masterlist=masterlist)
    

    # compact not isolated, nnsmall, lr
    M_small_nisol_S_nclustered_lr = M_small_nisol_S_nclustered.submask(lofarcat['LR'] >= lLR_thresh,
                        'lr',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), good LR'.format(s=size_large, nn=separation1),
                        edgelabel='Y',
                        qlabel='NN LR > {l:.2f}?'.format(l=lLR_thresh),
                        masterlist=masterlist)

    # compact not isolated, nnsmall, lr, NNlr
    M_small_nisol_S_nclustered_lr_NNlr = M_small_nisol_S_nclustered_lr.submask(lofarcat['NN_LR'] >= lLR_thresh,
                        'NNlr',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), good LR, NN good lr'.format(s=size_large, nn=separation1),
                        edgelabel='Y',
                        color='blue',
                        qlabel='accept LR',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_S_nclustered_lr_NNlr.mask] = 1

    # compact not isolated, nnsmall, lr, NNnlr
    M_small_nisol_S_nclustered_lr_NNnlr = M_small_nisol_S_nclustered_lr.submask(lofarcat['NN_LR'] < lLR_thresh,
                        'NNnlr',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), good LR, NN bad lr'.format(s=size_large, nn=separation1),
                        edgelabel='N',
                        color='blue',
                        qlabel='accept LR',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_S_nclustered_lr_NNnlr.mask] = 1



    # compact not isolated, nnsmall, nlr, NNnlr - there are possible doubles here!!
    M_small_nisol_S_nclustered_nlr_NNnlr = M_small_nisol_S_nclustered_nlr.submask(lofarcat['NN_LR'] < lLR_thresh,
                        'NNnlr',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), bad LR, NN bad lr'.format(s=size_large, nn=separation1),
                        edgelabel='N',
                        #color='orange',
                        qlabel='0.1 < flux ratio < 10?',
                        masterlist=masterlist)

    # compact not isolated, nnsmall, nlr, NNlr
    M_small_nisol_S_nclustered_nlr_NNlr = M_small_nisol_S_nclustered_nlr.submask(lofarcat['NN_LR'] >= lLR_thresh,
                        'NNlr',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), bad LR, NN good lr'.format(s=size_large, nn=separation1),
                        edgelabel='Y',
                        color='blue',
                        qlabel='accept no LR',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_S_nclustered_nlr_NNlr.mask] = 1

    C1_simflux = (lofarcat['NN_Frat'] <= 10) & (lofarcat['NN_Frat'] > 0.1)
    M_small_nisol_S_nclustered_nlr_NNnlr_simflux = M_small_nisol_S_nclustered_nlr_NNnlr.submask(C1_simflux,
                        'simflux',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), bad LR, NN bad lr, sim flux'.format(s=size_large, nn=separation1),
                        edgelabel='Y',
                        #color='blue',
                        qlabel='S1+S2 >= 50*(sep/100)**2 ?',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_S_nclustered_nlr_NNnlr_simflux.mask] = 1

    M_small_nisol_S_nclustered_nlr_NNnlr_diffflux = M_small_nisol_S_nclustered_nlr_NNnlr.submask(~C1_simflux,
                        'diffflux',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), bad LR, NN bad lr, diffflux'.format(s=size_large, nn=separation1),
                        edgelabel='N',
                        color='blue',
                        qlabel='accept no LR',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_S_nclustered_nlr_NNnlr_diffflux.mask] = 1

    C2_dist = ((lofarcat['NN_Total_flux']+lofarcat['Total_flux']) >= 50*(lofarcat['NN_sep']/100.)**2.)
    M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep = M_small_nisol_S_nclustered_nlr_NNnlr_simflux.submask(C2_dist,
                        'dist',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), bad LR, NN bad lr, sim flux'.format(s=size_large, nn=separation1),
                        edgelabel='Y',
                        qlabel='Bright?\n(S>{f2:.0f} mJy)'.format(f2=fluxcut2),
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep.mask] = 5
    
    
    
    M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_bright = M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep.submask(lofarcat['Total_flux'] >= fluxcut2,
                        'dist_bright',
                        'compact not isolated bright (s<{s:.0f}", NN<{nn:.0f}", S>= {f2:.0f}) NN small (s<={s:.0f}"), bad LR, NN bad lr, sim flux'.format(s=size_large, nn=separation1,f2=fluxcut2),
                        edgelabel='Y',
                        qlabel='ML output LGZ?',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_bright.mask] = 5
    
    
    M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_bright_mllgz = M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_bright.submask(lofarcat['ML_flag'] == False,
                        'mllgz',
                        'ml lgz',
                        edgelabel='Y',
                        qlabel='LGZ',
                        color='green',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_bright_mllgz.mask] = 3
    
    
    M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_bright_nmllgz = M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_bright.submask(lofarcat['ML_flag'] == True,
                        'nmllgz',
                        'not  ml lgz',
                        edgelabel='N',
                        qlabel='Visual sorting',
                        color='cyan',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_bright_nmllgz.mask] = 5
    
    
    
    
    M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_faint = M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep.submask(lofarcat['Total_flux'] < fluxcut2,
                        'dist_faint',
                        'compact not isolated faint (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}", S> {f2:.0f}), bad LR, NN bad lr, sim flux'.format(s=size_large, nn=separation1,f2=fluxcut2),
                        edgelabel='N',
                        qlabel='TBD',
                        color='orange',
                        masterlist=masterlist)
    lofarcat['ID_flag'][M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_faint.mask] = 5
    
    if step == 4:
        print ('not fully implemented')
        sys.exit(1)
        
        M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_lgz = M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep.submask(lofarcat['double_flag'] == 2,
                            'lgz',
                            edgelabel='complex',
                            color='green',
                            qlabel='lgz',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_lgz.mask] = 3
        
        M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_lgz3 = M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep.submask(lofarcat['double_flag2'] == 2,
                            'lgz3',
                            edgelabel='complex',
                            color='green',
                            qlabel='lgz3',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_lgz3.mask] = 330
        
        
        M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_noid = M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep.submask((lofarcat['double_flag'] == 1) | (lofarcat['double_flag2'] == 1),
                            'noid',
                            edgelabel='no match',
                            color='blue',
                            qlabel='accept no LR',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_noid.mask] = 1
        
        
        M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_art = M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep.submask((lofarcat['double_flag'] == 3) | (lofarcat['double_flag2'] == 3),
                            'artefact',
                            edgelabel='artefact',
                            color='grey',
                            qlabel='artefact',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_art.mask] = -1
        
        
        
        M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_prob = M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep.submask(((lofarcat['double_flag'] == 4) | (lofarcat['double_flag'] == 0)) & (lofarcat['double_flag2'] == 0),
                            'prob',
                            edgelabel='problem - noid',
                            color='blue',
                            qlabel='problem',
                            masterlist=masterlist)
        lofarcat['ID_flag'][M_small_nisol_S_nclustered_nlr_NNnlr_simflux_sep_prob.mask] = 1 ## these should be in lgz assoc! but let's double check...

    M_small_nisol_S_nclustered_nlr_NNnlr_simflux_nsep = M_small_nisol_S_nclustered_nlr_NNnlr_simflux.submask(~C2_dist,
                        'ndist',
                        'compact not isolated (s<{s:.0f}", NN<{nn:.0f}") NN small (s<={s:.0f}"), bad LR, NN bad lr, sim flux'.format(s=size_large, nn=separation1),
                        edgelabel='N',
                        color='blue',
                        qlabel='accept no LR',
                        masterlist=masterlist)    
    lofarcat['ID_flag'][M_small_nisol_S_nclustered_nlr_NNnlr_simflux_nsep.mask] = 1
    
    
    # other masks


    #M_encloses = Mask(lofarcat['Encloses'],
                    #'Encloses',
                    #'encloses',
                    #qlabel='encloses?')
    #M_enclosed = Mask(lofarcat['Enclosed'],
                    #'Enclosed',
                    #'enclosed',
                    #qlabel='enclosed?')
    #M_intersects = Mask(lofarcat['Intersects'],
                    #'Intersects',
                    #'intersects',
                    #qlabel='intersects?')
    #M_encloses.make_sample(lofarcat)
    #M_enclosed.make_sample(lofarcat)
    #M_intersects.make_sample(lofarcat)
    

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
    makesample = 0
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
        outname = 'flow_{h}_{version}_step{st:d}_s{s:.0f}_nn{nn:.0f}_flux{ff:.0f}'.format(h=h,version=version,st=step,s=size_large,nn=separation1,ff=fluxcut2)
        A.draw(outname+'.png') # write to file
        A.draw(outname+'.pdf') # write to file
        A.write(outname+'.dot') # write to file




    print('ID_flag count')
    t,i = np.unique(lofarcat['ID_flag'], return_counts=True)
    for tt,ii in zip(t,i): print(tt,ii)
    
    #print('ID_flag count - encloses')
    #t,i = np.unique(lofarcat['ID_flag'][lofarcat['Encloses']], return_counts=True)
    #for tt,ii in zip(t,i): print(tt,ii)
    #print('ID_flag count - enclosed')
    #t,i = np.unique(lofarcat['ID_flag'][lofarcat['Enclosed']], return_counts=True)
    #for tt,ii in zip(t,i): print(tt,ii)
    #print('ID_flag count - intersects')
    #t,i = np.unique(lofarcat['ID_flag'][lofarcat['Intersects']], return_counts=True)
    #for tt,ii in zip(t,i): print(tt,ii)
    
    
    print(fcflg,'count')
    t,i = np.unique(lofarcat[fcflg], return_counts=True)
    for tt,ii in zip(t,i): print(tt,ii)
    
    if 'MC_flag1' in lofarcat.colnames:
        print('MC_flag1 count')
        t,i = np.unique(lofarcat['MC_flag1'], return_counts=True)
        for tt,ii in zip(t,i): print(tt,ii)


    ## write output file
    # as hdf5 need to seralise the metadata
    lofarcat.write(lofarcat_file_srt, overwrite=True, serialize_meta=True)


    ## TESTING ##
    ## check gaus ML

