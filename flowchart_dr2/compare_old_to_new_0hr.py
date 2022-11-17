import numpy as np
from astropy.table import Table, join, vstack

#tpf = Table.read('/data2/wwilliams/projects/lofar_surveys/DR2/lgz_selection/LoTSS_DR2_{version}.srl_0h.lr-full.sorted_step2_flux4.prefilter_lgz_weave_selection_1.fits')
#length=6249
#told = Table.read('/data2/wwilliams/projects/lofar_surveys/DR2/LoTSS_DR2_{version}.srl_0h.lr-full.sorted_step2_flux4.fits')

#tcopy = Table.read('/data2/wwilliams/projects/lofar_surveys/DR2/flowchart_sep22/LoTSS_DR2_{version}.srl_0h.lr-full.sorted_step1_flux4.fits')


tsellgz = Table.read('/data2/wwilliams/projects/lofar_surveys/DR2/lgz_selection/LoTSS_DR2_{version}.srl_0h.lr-full.sorted_step2_flux4.lgz_selection_2.fits')
tselpf = Table.read('/data2/wwilliams/projects/lofar_surveys/DR2/lgz_selection/LoTSS_DR2_{version}.srl_0h.lr-full.sorted_step2_flux4.prefilter_lgz_selection_2.fits')

'''
t,i = np.unique(tsel_lgz1['FC_flag2'], return_counts=True)
for tt,ii in zip(t,i): print(tt,ii)
t,i = np.unique(tsel_lgz2['FC_flag2'], return_counts=True)
for tt,ii in zip(t,i): print(tt,ii)
t,i = np.unique(tsel_pf1['FC_flag2'], return_counts=True)
for tt,ii in zip(t,i): print(tt,ii)
t,i = np.unique(tsel_pf2['FC_flag2'], return_counts=True)
for tt,ii in zip(t,i): print(tt,ii)
'''
tsellgz_new = Table.read('/data2/wwilliams/projects/lofar_surveys/DR2/lgz_selection_nov18/LoTSS_DR2_{version}.srl_0h.lr-full.sorted_step2_flux4.lgz_selection_2.fits')
tselpf_new = Table.read('/data2/wwilliams/projects/lofar_surveys/DR2/lgz_selection_nov18/LoTSS_DR2_{version}.srl_0h.lr-full.sorted_step2_flux4.prefilter_lgz_selection_2.fits')

'''
t,i = np.unique(tsel_lgz1['FC_flag2'], return_counts=True)
for tt,ii in zip(t,i): print(tt,ii)
t,i = np.unique(tsel_lgz2['FC_flag2'], return_counts=True)
for tt,ii in zip(t,i): print(tt,ii)
t,i = np.unique(tsel_pf1['FC_flag2'], return_counts=True)
for tt,ii in zip(t,i): print(tt,ii)
t,i = np.unique(tsel_pf2['FC_flag2'], return_counts=True)
for tt,ii in zip(t,i): print(tt,ii)
'''

pfout = Table.read('LoTSS_DR2_{version}.srl_0h.prefilter_outputs.fits')
pfout = pfout[pfout['Prefilter']!=-99]      ## 10940

missing_lgz = []
lgz_selind = []
for s in tsellgz_new['Source_Name']:
    if s not in tsellgz['Source_Name']:
        missing_lgz.append(s)
        lgz_selind.append(True)
    else:
        lgz_selind.append(False)
not_needed_lgz = []
for s in tsellgz['Source_Name']:
    if s not in tsellgz_new['Source_Name']:
        not_needed_lgz.append(s)
        
missing_pf = []
for s in tselpf_new['Source_Name']:
    if s not in tselpf['Source_Name']:
        missing_pf.append(s)
missing_pf2 = []
pf_selind = []
for s in tselpf_new['Source_Name']:
    if s not in pfout['Source_Name']:
        missing_pf2.append(s)
        pf_selind.append(True)
    else:
        pf_selind.append(False)


tsellgz_newsel = tsellgz_new[lgz_selind]
tselpf_newsel = tselpf_new[pf_selind]

# remove any for lgz that already went to lgz through pf
lgz_selind = []
for s in tsellgz_newsel['Source_Name']:
    if s in pfout['Source_Name']:
        ii = (pfout['Source_Name'] == s)
        if pfout['Prefilter'][ii] in [1,4,7]:  # was sent to lgz / zoom / blend
            lgz_selind.append(False)
        else:
            lgz_selind.append(True)
    else:
        lgz_selind.append(True)
tsellgz_newsel_notpflgz = tsellgz_newsel[lgz_selind]

# remove any that have already gone to lgz:
pf_selind = []
for s in tselpf_newsel['Source_Name']:
    if s in tsellgz['Source_Name']:
        pf_selind.append(False)
    else:
        pf_selind.append(True)
tselpf_newsel_notlgz = tselpf_newsel[pf_selind]

print ('There are ',len(tsellgz_newsel),' sources selected this time for lgz that were not selected before for lgz')
print ('There are ',len(tsellgz_newsel_notpflgz),' sources selected this time for lgz that were not selected before for either lgz or pf (with lgz outcome)')
print ('There are ',len(tselpf_newsel),' sources selected this time for prefilter that were not selected before for prefilter')
print ('There are ',len(tselpf_newsel_notlgz),' sources selected this time for prefilter that were not selected before for either prefilter or lgz')


tsellgz_newsel_notpflgz.write('/data2/wwilliams/projects/lofar_surveys/DR2/lgz_selection_nov18/LoTSS_DR2_{version}.srl_0h.lr-full.sorted_step2_flux4.lgz_selection_2_fix.fits', overwrite=True)
tselpf_newsel.write('/data2/wwilliams/projects/lofar_surveys/DR2/lgz_selection_nov18/LoTSS_DR2_{version}.srl_0h.lr-full.sorted_step2_flux4.prefilter_lgz_selection_2_fix.fits', overwrite=True)

'''
lgz 1/2 old
5 3217
12 83
15 37
20 94
26 201
5 3172
12 79
15 54
20 75
26 295
lgz 1/2 new
5 3218
12 25
15 39
20 23
26 195
5 3172
12 14
15 54
20 12
26 287

prefilt 1/2 old
6 156
13 3075
16 31
21 2900
27 87
6 161
13 2739
16 15
21 1696
27 83

prefilt 1/2 new
6 156
13 698
16 31
21 1077
27 82
6 161
13 699
16 15
21 893
27 79
'''

