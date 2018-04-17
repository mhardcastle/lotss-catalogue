import os 
from astropy.table import Table, Column, join, vstack, hstack
from astropy.coordinates import SkyCoord
import numpy as np

from catalogue_create.process_lgz import Make_Shape

try:
    from tqdm import tqdm
    has_tqdm = True
except ImportError:
    has_tqdm = False

#from catalogue_create import process_lgz
'''
ID_flag
1 - ML
2 - 2MASX
3 - LGZ
32 - LGZ v2 pending
311 - lgz v1 - normal route
312 - lgz v1 zoom
4 - no id possible
5 - TBC
6 - deblending
61 - deblend directly
62 - deblend workflow
'''



'''
From the LOFAR catalogue:
Source_Name
RA
E_RA
E_RA_tot
DEC
E_DEC
E_DEC_tot
Peak_flux
E_Peak_flux
E_Peak_flux_tot
Total_flux
E_Total_flux
E_Total_flux_tot
Maj
E_Maj
Min
E_Min
PA
E_PA
Isl_rms
S_Code
Mosaic_ID
Isl_id

Append:
ID_flag
ID_wisename
ID_psname
ID_2masxname
ID_ra
ID_dec
ML_LR
LGZ_flags?

'''
def count_flags(cat, flag):
    print '{:s} counts'.format(flag)
    unique, counts = np.unique(cat[flag], return_counts=True)
    for u,c in zip(unique, counts):
        print u,c



def name_from_coords(ra,dec,prefix=''):
    sc = SkyCoord(ra,dec,frame='icrs', unit='deg')
    sc = sc.to_string(style='hmsdms',sep='',precision=2)
    name = prefix+sc.replace(' ','')[:-1]
    return name




if __name__=='__main__':

    ### Required INPUTS
    
    version =  '1.1'
    
    # lofar source catalogue
    path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'

    lofargcat_file = path+'LOFAR_HBA_T1_DR1_catalog_v0.99.gaus.fits'
    lofarcat_orig_file = path+'LOFAR_HBA_T1_DR1_catalog_v0.99.srl.gmasked.fits'

    # sorted output from flowchart
    lofarcat_file_srt = path+'LOFAR_HBA_T1_DR1_catalog_v0.99.srl.gmasked.sorted.fits'

    # LGZ output
    lgz_cat_file = os.path.join(path,'lgz_v2/HETDEX-LGZ-cat-v1.1-filtered-zooms.fits') 
    lgz_component_file = os.path.join(path,'lgz_v2/lgz_components.txt')

    comp_out_file = os.path.join(path,'LOFAR_HBA_T1_DR1_merge_ID_v{v:s}.comp.fits'.format(v=version))
    art_out_file = os.path.join(path,'LOFAR_HBA_T1_DR1_merge_ID_v{v:s}.art.fits'.format(v=version))
    merge_out_file = os.path.join(path,'LOFAR_HBA_T1_DR1_merge_ID_v{v:s}.fits'.format(v=version))
    merge_out_full_file = merge_out_file.replace('.fits','.full.fits')
    comp_out_full_file = comp_out_file.replace('.fits','.full.fits')

    srccat = Table.read(lofarcat_file_srt)
    lgz_components = Table.read(lgz_component_file, format='ascii', names=['Component_Name', 'Source_Name', 'lgz_flag'])
    
   
    lgz_cat_full = Table.read(lgz_cat_file)
    
    
    print 'Starting with {:d} sources'.format(len(srccat))

    ### add some needed columns
    srccat.add_column(Column(np.zeros(len(srccat),dtype='S60'),'ID_name'))
    srccat.add_column(Column(np.nan*np.zeros(len(srccat),dtype=float),'ID_ra'))
    srccat.add_column(Column(np.nan*np.zeros(len(srccat),dtype=float),'ID_dec'))
    
    ##
    srccat.add_column(Column(np.nan*np.zeros(len(srccat),dtype=float),'ML_LR'))

    srccat.add_column(Column(np.zeros(len(srccat),dtype=bool),'LGZ_art_flag'))
    
    

    # 2 LGZ sources are edge_flagged but not artefact in LGZ... manually set their LGZ Art_prob=1
    # so they will get removed from the component catalogue, and the src catalogue
    lgz_cat_full['Art_prob'][lgz_cat_full['Source_Name']=='ILTJ141331.77+530555.2'] = 1
    lgz_cat_full['Art_prob'][lgz_cat_full['Source_Name']=='ILTJ145944.21+455630.8'] = 1

    # 2 LGZ sources are 'deleted' in the zoom files...
    # do this manually rather than just dropping them at the end so they get in the artefact lists as well
    srccat['Artefact_flag'][srccat['Source_Name']=='ILTJ144521.43+472223.8'] = 1
    srccat['Artefact_flag'][srccat['Source_Name']=='ILTJ144953.51+483807.5'] = 1
    
    ## remove sources associated/flagged by LGZ
    # ideally this would just remove the components in the LGZ comp catalogue  - but using legacy catalogues mean that these don't directly map onto the new sources
    # martin has produced lgz_components.txt to do this.


    # make a temp table - rename CompName to SrcName to join to srccat
    tlgz_component = Table([lgz_components['Component_Name'],
                            Column(np.ones(len(lgz_components)),'LGZ_remove')])
    tlgz_component.rename_column('Component_Name','Source_Name')
    srccat.sort('Source_Name')
    tc = join(srccat, tlgz_component, join_type='left')
    tc['LGZ_remove'].fill_value = 0
    tc = tc.filled()
    tc.sort('Source_Name')
    srccat.add_column(tc['LGZ_remove'])
    del tc
    del tlgz_component


    is_lgz_art = np.zeros(len(srccat), dtype=bool)
    # keep a list of the lgz artefacts
    lgz_art = lgz_cat_full[lgz_cat_full['Art_prob'] >= 0.5]
    for ll in lgz_art:
        if ll['Assoc'] > 1:
            #get the comp naes
            llcomps = lgz_components[lgz_components['Source_Name'] ==  ll['Source_Name']]['Component_Name']
            for llcompi in llcomps:
                is_lgz_art[srccat['Source_Name'] == llcompi] = 1


        else:
            is_lgz_art[srccat['Source_Name'] == ll['Source_Name']] = 1

    print np.sum(is_lgz_art), 'lgz components are artefacts',
    print '(',len(lgz_art), 'lgz sources are artefacts)'
    
    srccat['LGZ_art_flag'][is_lgz_art] = 1
    srccat['Artefact_flag'][is_lgz_art] = 1

    # save list of the artefacts
    srccat[srccat['Artefact_flag'] == 1].write(art_out_file, overwrite=True)

    # keep the full component catalogue - renaming Source_Name to Component_Name
    compcat = srccat.copy()
    tc  = srccat['Source_Name'].copy()
    tc.name = 'Component_Name'
    compcat.add_column(tc)


    #### SRCCAT splits from COMPCAT here!! ###
    
    
    # keep track of the sources removed by LGZ - some of these need to be merged with 2MASX assoc
    lgz_select = (srccat['LGZ_remove']!=1)
    print 'Removing {n:d} sources associated in LGZ and LGZ artefacts'.format(n=np.sum(~lgz_select))
    srccat_del = srccat[~lgz_select]
    srccat = srccat[lgz_select]
    
    print 'getting the LGZ source names for the component table'
    # update the Source_Name's in the comp cat with the Source_Name's from the lgz_components:
    # first set them to LGZ  and there shouldn't be any of these left afterwards! will be check and crash if any are left
    # table merging is a pain but is much faster than numpy stuff
    compcat['Source_Name'][~lgz_select] = 'LGZ' 
    
    tcompcat = compcat.copy()
    tcompcat.rename_column('Source_Name','dummy') # need to match on Component_Name so get rid of the source_name
    tcompcat.sort('Component_Name')
    compcat.sort('Component_Name')
    tc = join(tcompcat, lgz_components, join_type='left', keys=['Component_Name'])
    tc.sort('Component_Name')
    assert len(tc) == len(compcat), 'table length mismatch after join'
    assert len(tc) == len(tcompcat), 'table length mismatch after join'
    for ii in [1,2]:
        selind = (tc['lgz_flag']==ii)
        compcat['Source_Name'][selind] = tc['Source_Name'][selind]
        compcat['ID_flag'][selind] = 310+ii
    del tc
    del tcompcat
        
    assert 'LGZ' not  in  compcat['Source_Name'], 'this should not happen! the LGZ sources have not all had their new names updated'
        
    ## remove artefacts
    # all the artefacts identified and visually confirmed in the flowchart process
    print 'Throwing away {n:d} artefacts'.format(n=np.sum(srccat['Artefact_flag'] == 1))
    srccat = srccat[srccat['Artefact_flag'] == 0]
    
    print 'Throwing away {n:d} component artefacts'.format(n=np.sum(compcat['Artefact_flag'] == 1))
    compcat = compcat[compcat['Artefact_flag'] == 0]
    

    # handle 2MASX sources
    ## HUGE 2MASX sources need to be removed, associated and added back
    ## the rest need a flag for 2MASX
    sel2mass = (srccat['ID_flag']==2)
    print 'adding info for {n:d} 2MASX source matches'.format(n=np.sum(sel2mass))
    # add the 2MASXJ
    names = srccat['2MASX_name'][sel2mass]
    names = np.array(['2MASX J'+n for n in names])
    
    
    srccat['ID_name'][sel2mass] = names
    srccat['ID_ra'][sel2mass] = srccat['2MASX_ra'][sel2mass]
    srccat['ID_dec'][sel2mass] = srccat['2MASX_dec'][sel2mass]
    
    
    # some sources come from a tag 'match to bright galaxy' - not necesarily 2MASX - look in SDSS for these:
    sel2masssdss = (srccat['ID_flag']==2) & (srccat['ID_name']=='2MASX J')
    
    print 'resorting to an SDSS match for {n:d} sources'.format(n=np.sum(sel2masssdss))
    
    import astropy.units as u
    from astroquery.sdss import SDSS
    snames = np.zeros(np.sum(sel2masssdss),dtype='S60')
    sdss_ra = np.zeros(np.sum(sel2masssdss))
    sdss_dec = np.zeros(np.sum(sel2masssdss))
    for ti,t in enumerate(srccat[sel2masssdss ]):
        ra,dec = t['RA'],t['DEC']
        #print ra,dec
        c = SkyCoord(ra,dec, frame='icrs', unit='deg')
        
        try:
            st = SDSS.query_region(c,radius=0.5*t['Maj']*u.arcsec, photoobj_fields=['ra','dec','objID','petroR50_r','petroMag_r'])
            st = st[(st['petroMag_r'] <20.) & (st['petroMag_r'] > 0 )] 
            #print st['petroMag_r'].max()
            c2 = SkyCoord(st['ra'],st['dec'], frame='icrs', unit='deg')
            sep = c.separation(c2)
            a = sep.argmin()
            sdss_ra[ti] = st['ra'][a]
            sdss_dec[ti] = st['dec'][a]
            snames[ti] = name_from_coords(sdss_ra[ti],sdss_dec[ti],prefix='SDSS J')
        except:
            print 'error - sdss', ra,dec
            import pdb ; pdb.set_trace()
            
        
    
    srccat['ID_name'][sel2masssdss] = snames
    srccat['ID_ra'][sel2masssdss] = sdss_ra
    srccat['ID_dec'][sel2masssdss] = sdss_dec
    
    srccat.add_column(Column(np.nan*np.zeros(len(srccat),dtype=float),'LGZ_Size'))
    srccat.add_column(Column(np.nan*np.zeros(len(srccat),dtype=float),'LGZ_Width'))
    srccat.add_column(Column(np.nan*np.zeros(len(srccat),dtype=float),'LGZ_PA'))
    srccat.add_column(Column(np.nan*np.zeros(len(srccat),dtype=float),'LGZ_Assoc'))
    
    
    srccat.add_column(Column(np.nan*np.zeros(len(srccat),dtype=float),'LGZ_Assoc_Qual'))
    srccat.add_column(Column(np.nan*np.zeros(len(srccat),dtype=float),'LGZ_ID_Qual'))

    print 'going to merge'
    
    ok_matches =  (names != '2MASXJ')
    
    unames, ucounts = np.unique(srccat['ID_name'][sel2mass ][ok_matches], return_counts=True)
    remove_2mass_mult = np.zeros(len(srccat),dtype=bool)
    lofarcat_add_2mass_mult = srccat[1:1]
    nmerge = np.sum(ucounts>1)
    nn = 0

    nok = 0
    nok1 = 0
    nnok = 0

    # we need to consider the 2MASX sources that were in the first LGZ selections (to be consistent with what came later we take the 2MASX merge and 'throw away' the lgz assoc as if they had never been in lgz) warning - this can be complicated if something is associated in lgz but not to the 2masx source... does this happen?
    # these are all the lgz sources that are mixed up with the 2msss sources - they should be removed from the lgz source list when these are added in later...
    # make sure all the components in the components cat point to the 2mass source
    remove_2mass_lgz_srcs = []
    
    if has_tqdm:
        iterate = tqdm(unames)
    else:
        iterate = unames
    for n in iterate: #[ucounts>1]:
        if 'SDSS' in n: continue
        mergeok = True
        
        i = np.where(srccat['2MASX_name'] == n.replace('2MASX J',''))[0]

        i2 = np.where(srccat_del['2MASX_name'] == n.replace('2MASX J',''))[0]
        
        if len(i2) > 0:
            nnok += 1
            srccat[i].write('check_large_opt/'+n.replace(' ','_')+'.fits',overwrite=True)

            lgz_srcs = []
            lgz_src_comps = []
            missing_comp = []
            for ii in i2:
                lgz_srci = lgz_components['Source_Name'][lgz_components['Component_Name'] == srccat_del['Source_Name'][ii]][0]
                all_lgz_comp = lgz_components['Component_Name'][lgz_components['Source_Name'] == lgz_srci]
                remove_2mass_lgz_srcs.append(lgz_srci)
                if lgz_srci not in lgz_srcs:
                    lgz_srcs.append(lgz_srci)
            for lgz_src in lgz_srcs:
                #print 'checking', lgz_src
                lgz_src_comp = lgz_components['Component_Name'][lgz_components['Source_Name'] == lgz_src]
                for lgz_src_compi in lgz_src_comp:
                    if lgz_src_compi not in lgz_src_comps:
                        lgz_src_comps.append(lgz_src_compi)

                    #if (lgz_src_compi not in srccat['Source_Name'][i])
                    if (lgz_src_compi not in srccat_del['Source_Name'][i2]):
                        if lgz_src_compi not in missing_comp:
                            missing_comp.append(lgz_src_compi)
                            print 'adding lgz component to 2MASX assoc', lgz_src_compi
                            i2 = np.hstack((i2, np.where(srccat_del['Source_Name'] == lgz_src_compi)[0] ))
                            mergeok = False

            srccat_del[i2].write('check_large_opt/'+n.replace(' ','_')+'-lgz.fits',overwrite=True)
        else:
            if len(i) == 1:
                nok1 += 1
            else:
                nok += 1

        nn += len(i)
        remove_2mass_mult[i] = True
        
        complist = vstack(( srccat[i].copy() , srccat_del[i2].copy() ))
        assoc_2mass = srccat[i].copy()[0]
        
        
        # merging multiple S/M will be M, unless merging 1 S source
        if len(complist) > 0:
            
            assoc_2mass['RA']=np.average(complist['RA'], weights=complist['Total_flux'])
            assoc_2mass['DEC']=np.average(complist['DEC'], weights=complist['Total_flux'])
            
            assoc_2mass['ID_flag'] = 2 # ensure the ID_flag is 2 (if the first of the other components was not 2)
            #  (if the first of the other components was not 2 then it won't have the right ID_name)
            if np.any(complist['ID_name'] == ''):
                assert len(np.unique(complist['ID_name'][complist['ID_name']!=''])), 'something is wrong, there are multiple 2MASX ids'
                assoc_2mass['ID_name'] = complist['ID_name'][complist['ID_name']!=''][0]
                
                
            
            assoc_2mass['Source_Name'] = name_from_coords(assoc_2mass['RA'],assoc_2mass['DEC'],prefix='ILTJ')
            
            assoc_2mass['E_RA']=np.sqrt(np.sum(complist['E_RA']**2.0))/len(complist)
            assoc_2mass['E_DEC']=np.sqrt(np.sum(complist['E_DEC']**2.0))/len(complist)
            assoc_2mass['Isl_rms']=np.mean(complist['Isl_rms'])
            assoc_2mass['Total_flux']=np.sum(complist['Total_flux'])
            # total flux error is error on the sum
            assoc_2mass['E_Total_flux']=np.sqrt(np.sum(complist['E_Total_flux']**2.0))
            # peak flux and error from brightest component
            maxpk=np.argmax(complist['Peak_flux'])
            assoc_2mass['Peak_flux']=complist[maxpk]['Peak_flux']
            assoc_2mass['E_Peak_flux']=complist[maxpk]['E_Peak_flux']
            
            assoc_2mass['Isl_id'] = -99
            if len(complist) > 1:
                assoc_2mass['S_Code'] = 'M'
            #else:
                #assoc_2mass['S_Code'] = 'S'
            
        for t in ['Maj', 'Min', 'PA']:
            assoc_2mass[t] = np.nan
            assoc_2mass['E_'+t] = np.nan
            assoc_2mass['DC_'+t] = np.nan
            assoc_2mass['E_DC_'+t] = np.nan
                        
        
        #c =SkyCoord(complist['RA'], complist['DEC'], unit='deg')
        # TBD 'Mosiac_ID'
        
        # size is taken from convex hull of components - as in LGZ process
        if len(complist) > 0:
            cshape = Make_Shape(complist)
            assoc_2mass['LGZ_Size'] = cshape.length()
            assoc_2mass['LGZ_Width'] = cshape.width()
            assoc_2mass['LGZ_PA'] = cshape.pa()
            
            #if cshape.length() < 1:
                #print complist
            #print np.array(complist['Source_Name']), np.array(complist['Maj']), cshape.length()
        else:
            assoc_2mass['LGZ_Size'] = complist['DC_Maj'][0]
            assoc_2mass['LGZ_Width'] = complist['DC_Min'][0]
            assoc_2mass['LGZ_PA'] = complist['DC_PA'][0]
        
        
        assoc_2mass['LGZ_Assoc'] = len(complist)
        
        # give them quality flags like LGZ
        assoc_2mass['LGZ_Assoc_Qual'] = 1.
        assoc_2mass['LGZ_ID_Qual'] = 1.
                
        # to save the new names
        for c in complist:
            ni = np.where(compcat['Component_Name'] == c['Source_Name'])[0]
            compcat['Source_Name'][ni] = assoc_2mass['Source_Name']
            compcat['ID_flag'][ni] = 2
        
        if assoc_2mass['ID_name'] == '':
            sys.exit()
        
        lofarcat_add_2mass_mult = vstack([lofarcat_add_2mass_mult, assoc_2mass]) 
        
    print nok, 'are ok'
    print nok1, 'are ok (single)'
    print nnok, 'are partly in an lgz source'
    #sys.exit()
    srccat = srccat[~remove_2mass_mult]
    srccat = vstack([srccat,lofarcat_add_2mass_mult])
        
    print 'merging components of {n:d} 2MASX sources'.format(n=nmerge)
    print 'removing merged components {n:d} from catalogue'.format(n=np.sum(remove_2mass_mult))
    print 'adding back {n:d} merged sources'.format(n=len(lofarcat_add_2mass_mult))
    
    
    
    # handle ML sources
    lLR_thresh = 0.639
    selml = ((srccat['ID_flag']==1) |(srccat['ID_flag']==61) | (srccat['ID_flag']==62)) & (np.log10(1+srccat['LR']) > lLR_thresh)
    print 'adding info for {n:d} ML source matches'.format(n=np.sum(selml))
    

    # take the PS name over the WISE name
    # why is PS name just some number ?? - pepe?
    namesP = srccat['LR_name_ps'][selml]
    namesP = [ 'PS '+str(nP) if nP != 999999 else '' for nP in namesP ]
    #namesP = [name_from_coords(ra,dec, prefix='PSO J') for ra,dec,n in zip(srccat['LR_ra'][selml],srccat['LR_dec'][selml],srccat['LR_name_ps'][selml])  ]
    namesW = srccat['LR_name_wise'][selml]
    namesW = [ 'AllWISE'+nW  if nW != 'N/A' else '' for nW in namesW]
    names = [nP if nP != '' else nW  for nP,nW in zip(namesP,namesW)]
    
    
    srccat['ID_name'][selml] = names
    srccat['ID_ra'][selml] = srccat['LR_ra'][selml]
    srccat['ID_dec'][selml] = srccat['LR_dec'][selml]
    srccat['ML_LR'][selml] = srccat['LR'][selml]
    
    # use gaus info where it is needed:
    # for the blends where there is no source match
    selmlg_blend = ((srccat['ID_flag']==61) | (srccat['ID_flag']==62)) & (np.log10(1+srccat['LR']) <= lLR_thresh)
    # and for the msources that are selected to have best G match
    selmlg_auto = (srccat['ID_flag']==1) & ( ((srccat['msource1_flag']==2) | (srccat['msource2_flag']==2)))
    selmlg = selmlg_blend | selmlg_auto
    print 'adding info for {n:d} ML gaus source matches'.format(n=np.sum(selmlg))
    

    # take the PS name over the WISE name
    # why is PS name just some number ?? - pepe?
    namesP = srccat['gLR_name_ps'][selmlg]
    namesP = [ 'PS '+str(nP) if nP != 999999 else '' for nP in namesP ]
    #namesP = [name_from_coords(ra,dec, prefix='PSO J') for ra,dec,n in zip(srccat['LR_ra'][selmlg],srccat['LR_dec'][selmlg],srccat['LR_name_ps'][selmlg])  ]
    namesW = srccat['gLR_name_wise'][selmlg]
    namesW = [ 'AllWISE'+nW  if nW != 'N/A' else '' for nW in namesW]
    names = [nP if nP != '' else nW  for nP,nW in zip(namesP,namesW)]
    
    
    srccat['ID_name'][selmlg] = names
    srccat['ID_ra'][selmlg] = srccat['gLR_ra'][selmlg]
    srccat['ID_dec'][selmlg] = srccat['gLR_dec'][selmlg]
    srccat['ML_LR'][selmlg] = srccat['gLR'][selmlg]
    
    
    
    selml = (srccat['ID_flag']==1) & (np.log10(1+srccat['LR']) <= lLR_thresh)
    print 'adding info for {n:d} ML source non-matches'.format(n=np.sum(selml))
    
    srccat['ID_name'][selml] = ''

    lgz_cat_full.add_column(Column(np.zeros(len(lgz_cat_full), dtype=bool),'2MASSoverlap'))
    for s in remove_2mass_lgz_srcs:
        lgz_cat_full['2MASSoverlap'][lgz_cat_full['Source_Name'] == s] = True
    print np.sum(lgz_cat_full['2MASSoverlap']), ' removed because they are in a 2MASS assoc'

    # get a list of the lgz artefacts...
    lgz_art['Assoc'] > 1

    ## add LGz v1 associated sources
    # 
    lgz_select = (lgz_cat_full['Compoverlap']==0)&(lgz_cat_full['Art_prob']<0.5)&(lgz_cat_full['Zoom_prob']<0.5)&(lgz_cat_full['2MASSoverlap']==0)
    print 'Selecting {n2:d} of {n1:d} sources in the LGZ catalogue to add'.format(n1=len(lgz_cat_full),n2=np.sum(lgz_select))
    lgz_cat = lgz_cat_full[lgz_select]
    lgz_cat.rename_column('optRA','ID_ra')
    lgz_cat.rename_column('optDec','ID_dec')
    lgz_cat.rename_column('OptID_Name','ID_name')
    #lgz_cat.rename_column('Size','LGZ_Size')
    lgz_cat.rename_column('Assoc','LGZ_Assoc')
    lgz_cat.rename_column('Assoc_Qual','LGZ_Assoc_Qual')
    lgz_cat.rename_column('ID_Qual','LGZ_ID_Qual')
    
    lgz_cat.rename_column('New_size','LGZ_Size')
    lgz_cat.rename_column('New_width','LGZ_Width')
    lgz_cat.rename_column('New_PA','LGZ_PA')
    
    lgz_cat.add_column(Column(3*np.ones(len(lgz_cat),dtype=int),'ID_flag'))
    for lgzci in lgz_components:
        ls, lf  =  lgzci['Source_Name'], lgzci['lgz_flag']
        ind = np.where(lgz_cat['Source_Name'] == ls)[0]
        if lf == 1:
            lgz_cat['ID_flag'][ind] = 311
        elif lf == 2:
            lgz_cat['ID_flag'][ind] = 312
        else:
            print 'error'
            
        
    
    ## change None to ''
    lgz_cat['ID_name'][lgz_cat['ID_name']=='None'] = ''

    ## get rid of the strange metadata in the lgz cat...
    lgz_cat.meta = srccat.meta
    mergecat = vstack([srccat, lgz_cat])
    print 'now we have {n:d} sources'.format(n=len(mergecat))

    
    
    ### this is blend and should be dealt with later:
    # but update the ID_flag before we discard the Blend_prob column
    mergecat['ID_flag'][mergecat['Blend_prob'] >= 0.5] = 63
    
    # update ID_flags of the blend components from the component file
    lgz_select_blend = (lgz_cat_full['Blend_prob']>=0.5)
    print 'Selecting {n2:d} of {n1:d} blend sources in the LGZ catalogue to update ID_flag'.format(n1=len(lgz_cat_full),n2=np.sum(lgz_select_blend))
    lgz_cat_blend = lgz_cat_full[lgz_select_blend]
    for s in lgz_cat_blend['Source_Name']:
        #if s not in mergecat['Source_Name'][mergecat['Blend_prob'] >= 0.5]:
            #if s in mergecat['Source_Name']:
                #print s, mergecat['ID_flag','Zoom_prob','Blend_prob'][mergecat['Source_Name']==s][0]
        # look up component names
        si = (lgz_components['Source_Name'] == s)
        lgz_c = lgz_components['Component_Name'][si]
        for c in lgz_c:
            ci = compcat['Component_Name'] == c
            compcat['ID_flag'][ci] = 63
    
    
    
    # write some flag counts for both catalogues    
    count_flags(mergecat, 'ID_flag')
    count_flags(compcat, 'ID_flag')
    
    compcat.write(comp_out_full_file, overwrite=True)
    
    compcat.keep_columns(['Component_Name', 'Source_Name', 'RA', 'E_RA', 'DEC', 'E_DEC', 'Peak_flux', 'E_Peak_flux', 'Total_flux', 'E_Total_flux', 'Maj', 'E_Maj', 'Min', 'E_Min', 'PA', 'E_PA', 'DC_Maj', 'E_DC_Maj', 'DC_Min', 'E_DC_Min', 'DC_PA', 'E_DC_PA', 'Isl_rms', 'S_Code', 'Ng', 'Mosaic_ID', 'Number_Masked', 'Number_Pointings', 'Masked_Fraction', 'ID_flag'])
    
    compcat.write(comp_out_file, overwrite=True)

    mergecat.write(merge_out_full_file, overwrite=True)


    ## throw away extra columns
    mergecat.keep_columns(['Source_Name', 'RA', 'E_RA', 'DEC', 'E_DEC', 'Peak_flux', 'E_Peak_flux', 'Total_flux', 'E_Total_flux', 'Maj', 'E_Maj', 'Min', 'E_Min', 'PA', 'E_PA', 'DC_Maj', 'E_DC_Maj', 'DC_Min', 'E_DC_Min', 'DC_PA', 'E_DC_PA', 'Isl_rms', 'S_Code', 'Mosaic_ID', 'Number_Masked', 'Number_Pointings', 'Masked_Fraction', 'ID_flag', 'ID_name', 'ID_ra', 'ID_dec', 'ML_LR', 'LGZ_Size', 'LGZ_Width', 'LGZ_PA', 'LGZ_Assoc', 'LGZ_Assoc_Qual', 'LGZ_ID_Qual'])

    mergecat.write(merge_out_file, overwrite=True)
    