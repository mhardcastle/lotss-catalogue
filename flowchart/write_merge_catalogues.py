import os 
from astropy.table import Table, Column, join, vstack
from astropy.coordinates import SkyCoord
import numpy as np
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






if __name__=='__main__':

    ### Required INPUTS
    # lofar source catalogue, gaussian catalogue and ML catalogues for each


    path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'

    lofargcat_file = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.gaus.fixed.fits'
    lofarcat_orig_file = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.srl.fixed.fits'

    # PS ML - matches for sources and gaussians
    psmlcat_file = path+'lofar_pw.fixed.fits'
    psmlgcat_file = path+'lofar_gaus_pw.fixed.fits'

    # sorted output from flowchart
    lofarcat_file_srt = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.srl.fixed.sorted.fits'

    # LGZ output
    #lgz_compcat_file = os.path.join(path,'LGZ_v0/HETDEX-LGZ-comps-v0.5.fits')
    #lgz_cat_file = os.path.join(path,'LGZ_v0/HETDEX-LGZ-cat-v0.5-filtered.fits') 
    lgz_cat_file = os.path.join(path,'lgz_v1/HETDEX-LGZ-cat-v0.6-filtered-zooms.fits') 
    lgz_component_file = os.path.join(path,'lgz_v1/lgz_components.txt')

    comp_out_file = os.path.join(path,'LOFAR_HBA_T1_DR1_merge_ID_v0.6.comp.fits')
    merge_out_file = os.path.join(path,'LOFAR_HBA_T1_DR1_merge_ID_v0.6.fits')
    merge_out_full_file = merge_out_file.replace('.fits','.full.fits')

    lofarcat_sorted = Table.read(lofarcat_file_srt)
    lofarcat_sorted_antd = Table.read(lofarcat_file_srt)
    with open(lgz_component_file,'r') as f:
        lgz_lines = f.readlines()
    lgz_component = [l.rstrip().split()[0] for l in lgz_lines]
    lgz_src = [l.rstrip().split()[1] for l in lgz_lines]
    lgz_flag = [int(l.rstrip().split()[2]) for l in lgz_lines]
    
    #psmlcat = Table.read(psmlcat_file)
    
    #lgz_compcat = Table.read(lgz_compcat_file)
    lgz_cat = Table.read(lgz_cat_file)
    
    
    print 'Starting with {:d} sources'.format(len(lofarcat_sorted))

    ### add some needed columns
    
    lofarcat_sorted.add_column(Column(np.zeros(len(lofarcat_sorted),dtype='S60'),'ID_name'))
    #lofarcat_sorted.add_column(Column(['None']*len(lofarcat_sorted),'ID_name'))
    lofarcat_sorted.add_column(Column(np.nan*np.zeros(len(lofarcat_sorted),dtype=float),'ID_ra'))
    lofarcat_sorted.add_column(Column(np.nan*np.zeros(len(lofarcat_sorted),dtype=float),'ID_dec'))
    
    
    ##
    lofarcat_sorted.add_column(Column(np.nan*np.zeros(len(lofarcat_sorted),dtype=float),'ML_LR'))

    #lofarcat_sorted_antd.add_column(Column(np.nan*np.zeros(len(lofarcat_sorted_antd),dtype=float),'Removed'))
    tc  = lofarcat_sorted['Source_Name'].copy()
    tc.name = 'New_Source_Name'
    lofarcat_sorted_antd.add_column(tc)

    ## remove sources associated/flagged by LGZ v0
    # ideally this would just remove the components in the LGZ comp catalogue  - but using legacy catalogues mean that these don't directly map onto the new sources
    # martin has produced remove.txt to do this.


    lgz_select = np.ones(len(lofarcat_sorted), dtype=bool)
    #for si,s in enumerate(lofarcat_sorted['Source_Name']):
        #if s in lgz_component:
            #lgz_select[si] = False
    lgz_component = np.unique(lgz_component)
    tlgz_component = Table([Column(lgz_component,'Source_Name'), Column(np.ones(len(lgz_component)),'LGZ_remove')])
    lofarcat_sorted.sort('Source_Name')
    tc = join(lofarcat_sorted, tlgz_component, join_type='left')
    tc['LGZ_remove'].fill_value = 0
    tc = tc.filled()
    tc.sort('Source_Name')
    lgz_select = (tc['LGZ_remove']!=1)

    #import ipdb ; ipdb.set_trace()

    print 'Removing {n:d} sources associated in LGZ v1'.format(n=np.sum(~lgz_select))
    lofarcat_sorted = lofarcat_sorted[lgz_select]
    # we don't know what their new names are
    lofarcat_sorted_antd['New_Source_Name'][~lgz_select] = 'LGZ' # there shouldn't be any of these left afterwards!
    for lc, ls, lf in zip(lgz_component, lgz_src, lgz_flag):
        ind = np.where(lofarcat_sorted_antd['Source_Name'] == lc)[0]
        lofarcat_sorted_antd['New_Source_Name'][ind] = ls
        
        if lf == 1:
            lofarcat_sorted_antd['ID_flag'][ind] = 311
        elif lf == 2:
            lofarcat_sorted_antd['ID_flag'][ind] = 312

    ## remove artefacts
    # all the artefacts identified and visually confirmed in the flowchart process
    print 'Throwing away {n:d} artefacts'.format(n=np.sum(lofarcat_sorted['Artefact_flag'] == 1))
    lofarcat_sorted = lofarcat_sorted[lofarcat_sorted['Artefact_flag'] == 0]
    
    # artefacts have no name in the merged catalogue cos they don't appear there
    lofarcat_sorted_antd['New_Source_Name'][lofarcat_sorted_antd['Artefact_flag'] != 0] = ''
    
    print 'left with {n:d} sources'.format(n=len(lofarcat_sorted))




    # handle TBC

    # handle 2MASX sources
    ## HUGE 2MASX sources need to be removed, associated and added back
    ## the rest need a flag for 2MASX
    sel2mass = (lofarcat_sorted['ID_flag']==2)
    print 'adding info for {n:d} 2MASX source matches'.format(n=np.sum(sel2mass))
    # add the 2MASXJ
    names = lofarcat_sorted['2MASX_name'][sel2mass]
    names = np.array(['2MASXJ'+n for n in names])
    
    
    lofarcat_sorted['ID_name'][sel2mass] = names
    lofarcat_sorted['ID_ra'][sel2mass] = lofarcat_sorted['2MASX_ra'][sel2mass]
    lofarcat_sorted['ID_dec'][sel2mass] = lofarcat_sorted['2MASX_dec'][sel2mass]
    
    
    # some sources come from a tag 'match to bright galaxy' - not necesarily 2MASX - look in SDSS for these:
    sel2masssdss = (lofarcat_sorted['ID_flag']==2) & (lofarcat_sorted['ID_name']=='2MASXJ')
    #sdss_matches =  (names == '2MASXJ')
    
    #lofarcat_sorted['ID_name'][sel2mass ] = names
    #lofarcat_sorted[sel2mass ][sdss_matches]
    
    print 'resorting to an SDSS match for {n:d} sources'.format(n=np.sum(sel2masssdss))
    
    #### TBD ####
    import astropy.units as u
    from astroquery.sdss import SDSS
    snames = np.zeros(np.sum(sel2masssdss),dtype='S60')
    sdss_ra = np.zeros(np.sum(sel2masssdss))
    sdss_dec = np.zeros(np.sum(sel2masssdss))
    for ti,t in enumerate(lofarcat_sorted[sel2masssdss ]):
        ra,dec = t['RA'],t['DEC']
        #print ra,dec
        c = SkyCoord(ra,dec, frame='icrs', unit='deg')
        
        #c.
        
        try:
            st = SDSS.query_region(c,radius=0.5*t['Maj']*u.arcsec, photoobj_fields=['ra','dec','objID','petroR50_r','petroMag_r'])
            st = st[(st['petroMag_r'] <20.) & (st['petroMag_r'] > 0 )] 
            #print st['petroMag_r'].max()
            c2 = SkyCoord(st['ra'],st['dec'], frame='icrs', unit='deg')
            sep = c.separation(c2)
            a = sep.argmin()
            #print st
            #print st[a]
            snames[ti] = 'SDSS '+str(st['objID'][a])
            sdss_ra[ti] = st['ra'][a]
            sdss_dec[ti] = st['dec'][a]
        except:
            print 'error - sdss', ra,dec
        
    
    #c = SkyCoord(lofarcat_sorted['RA'][sel2mass ][sdss_matches], lofarcat_sorted['DEC'][sel2mass ][sdss_matches], frame='icrs', unit='deg')
    #sdss_tab = SDSS.query_crossid(c,radius=20*u.arcsec)
    
    lofarcat_sorted['ID_name'][sel2masssdss] = snames
    lofarcat_sorted['ID_ra'][sel2masssdss] = sdss_ra
    lofarcat_sorted['ID_dec'][sel2masssdss] = sdss_dec
    
    lofarcat_sorted.add_column(Column(np.nan*np.zeros(len(lofarcat_sorted),dtype=float),'LGZ_Size'))
    lofarcat_sorted.add_column(Column(np.nan*np.zeros(len(lofarcat_sorted),dtype=float),'LGZ_Assoc'))

    ok_matches =  (names != '2MASXJ')
    
    unames, ucounts = np.unique(lofarcat_sorted['ID_name'][sel2mass ][ok_matches], return_counts=True)
    remove_2mass_mult = np.zeros(len(lofarcat_sorted),dtype=bool)
    lofarcat_add_2mass_mult = lofarcat_sorted[1:1]
    nmerge = np.sum(ucounts>1)
    nn = 0
    for n in unames[ucounts>1]:
        
        i = np.where(lofarcat_sorted['ID_name'] == n)[0]
        nn += len(i)
        #print n, i
        remove_2mass_mult[i] = True
        
        complist = lofarcat_sorted[i]
        assoc_2mass = lofarcat_sorted[i[0]]
        
        assoc_2mass['RA']=np.average(complist['RA'], weights=complist['Total_flux'])
        assoc_2mass['DEC']=np.average(complist['DEC'], weights=complist['Total_flux'])
        sc = SkyCoord(assoc_2mass['RA'],assoc_2mass['DEC'],frame='icrs', unit='deg')
        sc = sc.to_string(style='hmsdms',sep='',precision=2)
        assoc_2mass['Source_Name'] = str('ILTJ'+sc).replace(' ','')[:-1]
        
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
        # merging multiple S/M will be M
        assoc_2mass['S_Code'] = 'M'
        for t in ['Maj', 'Min', 'PA','E_Maj', 'E_Min', 'E_PA']:
            assoc_2mass[t] = np.nan
        assoc_2mass['Isl_id'] = -99
        
        c =SkyCoord(complist['RA'], complist['DEC'], unit='deg')
                

        # size is max of Maj or sep between centres
        assoc_2mass['LGZ_Size']=np.max([np.max(complist['Maj']) , np.max([ci.separation(c).max().to('arcsec').value for ci in c])])
        # TBD 'Mosiac_ID'
        assoc_2mass['LGZ_Assoc'] = len(complist)
        
        #import ipdb ; ipdb.set_trace()
        
        # to save the new names
        for c in complist:
            ni = np.where(lofarcat_sorted_antd['Source_Name'] == c['Source_Name'])[0]
            lofarcat_sorted_antd['New_Source_Name'][ni] = assoc_2mass['Source_Name']
        
        
        lofarcat_add_2mass_mult = vstack([lofarcat_add_2mass_mult, assoc_2mass]) 
        
    lofarcat_sorted = lofarcat_sorted[~remove_2mass_mult]
    lofarcat_sorted = vstack([lofarcat_sorted,lofarcat_add_2mass_mult])
        
    print 'merging components of {n:d} 2MASX sources'.format(n=nmerge)
    print 'removing merged components {n:d} from catalogue'.format(n=np.sum(remove_2mass_mult))
    print 'adding back {n:d} merged sources'.format(n=len(lofarcat_add_2mass_mult))
    
    
    
    # handle ML sources
    lLR_thresh = 0.36
    selml = (lofarcat_sorted['ID_flag']==1) & (np.log10(1+lofarcat_sorted['LR']) > lLR_thresh)
    print 'adding info for {n:d} ML source matches'.format(n=np.sum(selml))
    

    
    # take the PS name over the WISE name
    # why is PS name just some number ?? - pepe?
    namesP = lofarcat_sorted['LR_name_ps'][selml]
    namesW = lofarcat_sorted['LR_name_wise'][selml]
    names = [ 'PS '+str(nP) if nP != 999999  else 'AllWISE'+nW for nP,nW in zip(namesP,namesW)]
    
    
    lofarcat_sorted['ID_name'][selml] = names
    lofarcat_sorted['ID_ra'][selml] = lofarcat_sorted['LR_ra'][selml]
    lofarcat_sorted['ID_dec'][selml] = lofarcat_sorted['LR_dec'][selml]
    lofarcat_sorted['ML_LR'][selml] = lofarcat_sorted['LR'][selml]
    
    selml = (lofarcat_sorted['ID_flag']==1) & (np.log10(1+lofarcat_sorted['LR']) <= lLR_thresh)
    print 'adding info for {n:d} ML source non-matches'.format(n=np.sum(selml))
    
    lofarcat_sorted['ID_name'][selml] = ''

                               

    ## add LGz v0 associated sources
    # 
    lgz_select = (lgz_cat['Compoverlap']==0)&(lgz_cat['Art_prob']<0.5)&(lgz_cat['Zoom_prob']<0.5)&(lgz_cat['Blend_prob']<0.5)&(lgz_cat['Hostbroken_prob']<0.5)
    print 'Selecting {n2:d} of {n1:d} sources in the LGZv0 catalogue to add'.format(n1=len(lgz_cat),n2=np.sum(lgz_select))
    lgz_cat = lgz_cat[lgz_select]
    lgz_cat.rename_column('optRA','ID_ra')
    lgz_cat.rename_column('optDec','ID_dec')
    lgz_cat.rename_column('OptID_Name','ID_name')
    lgz_cat.rename_column('Size','LGZ_Size')
    lgz_cat.rename_column('Assoc','LGZ_Assoc')
    lgz_cat.rename_column('Assoc_Qual','LGZ_Assoc_Qual')
    lgz_cat.rename_column('ID_Qual','LGZ_ID_Qual')
    
    lgz_cat.add_column(Column(3*np.ones(len(lgz_cat),dtype=int),'ID_flag'))
    for ls, lf in zip(lgz_src, lgz_flag):
        ind = np.where(lgz_cat['Source_Name'] == ls)[0]
        if lf == 1:
            lgz_cat['ID_flag'][ind] = 311
        elif lf == 2:
            lgz_cat['ID_flag'][ind] = 312
        else:
            print 'error'
        

    ## change None to ''
    lgz_cat['ID_name'][lgz_cat['ID_name']=='None'] = ''

    mergecat = vstack([lofarcat_sorted, lgz_cat])
    print 'now we have {n:d} sources'.format(n=len(mergecat))

    print 'ID_flag counts:'
    unique, counts = np.unique(mergecat['ID_flag'], return_counts=True)
    for u,c in zip(unique, counts):
        print u,c


    if os.path.isfile(comp_out_file):
        os.remove(comp_out_file)
    lofarcat_sorted_antd.write(comp_out_file)

    if os.path.isfile(merge_out_full_file):
        os.remove(merge_out_full_file)
    mergecat.write(merge_out_full_file)


    ## throw away extra columns
    mergecat.keep_columns(['Source_Name', 'RA', 'E_RA', 'DEC', 'E_DEC', 'Peak_flux', 'E_Peak_flux', 'Total_flux', 'E_Total_flux', 'Maj', 'E_Maj', 'Min', 'E_Min', 'PA', 'E_PA', 'Isl_rms', 'S_Code', 'Mosaic_ID', 'ID_flag', 'ID_name', 'ID_ra', 'ID_dec', 'ML_LR', 'LGZ_Size', 'LGZ_Assoc', 'LGZ_Assoc_Qual', 'LGZ_ID_Qual'])

    
    if os.path.isfile(merge_out_file):
        os.remove(merge_out_file)
    mergecat.write(merge_out_file)
    
    sys.exit()
    tt = mergecat['ID_name']
    tt = tt[tt!='']
    tt = tt[tt!='2MASXJ']
    tt = tt[tt!='Mult']
    n,u=np.unique(tt,return_counts=True)

    tout = mergecat[1:1]
    for nn,uu in zip(n,u):
        if uu > 1: 
            #print nn,uu
            tt=mergecat[mergecat['ID_name']==nn]
            if tt['ID_flag'][0] == 1:
                xx  = Table([tt['Source_Name'], tt['ID_flag'], tt['ID_name'], tt['FC_flag'], tt['Art_prob'],tt['Blend_prob'],tt['Zoom_prob'],tt['Hostbroken_prob'],tt['Total_flux']])
                print xx
                tout = vstack([tout,tt])
