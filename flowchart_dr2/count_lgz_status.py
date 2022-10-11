import sys
import numpy as np
from astropy.table import Table


import locale
locale.setlocale(locale.LC_ALL, 'en_US')  # for , separators


#print(len(cat),'sources in catalogue',catname)
#print(fcflg,'count')
#m = cat['WEAVE_priority{}'.format(weave_pri)]==True


size_large = 15.           # in arcsec
separation1 = 45.          # in arcsec
size_huge = 25.            # in arcsec
#separation2 = 30.          # in arcsec
fluxcut = 8.               # in mJy
fluxcut2 = 4.   #float(sys.argv[3])


step=2
path = '/Users/w.williams/projects/lofar_surveys/DR2/'
version ='v110'

cat0 = Table.read(path+'LoTSS_DR2_{version}.srl_{h}.lr-full.sorted_step{st}_flux{ff:.0f}.hdf5'.format(h='0h',version=version,st=step,ff=fluxcut2))
cat13 = Table.read(path+'LoTSS_DR2_{version}.srl_{h}.lr-full.sorted_step{st}_flux{ff:.0f}.hdf5'.format(h='13h',version=version,st=step,ff=fluxcut2))
#cat13d = Table.read(path+'LoTSS_DR2_{version}.srl_{h}.lr-full.sorted_step{st}_flux{ff:.0f}_weavepri12.hdf5'.format(h='13h',version=version,st=step,ff=fluxcut2))

N0 = len(cat0)
N13 = len(cat13)

print('all')
print ('0hr: {0:d}'.format(N0))
print ('13hr: {0:d}'.format(N13))

print()
print ('large & bright')
N0 = np.sum((cat0['Total_flux']>=fluxcut)&(cat0['Maj']>size_large))
N13 = np.sum((cat13['Total_flux']>=fluxcut)&(cat13['Maj']>size_large))
print ('0hr: {0:d}'.format(N0))
print ('13hr: {0:d}'.format(N13))

print()
print('WEAVE_priority 1')
N0 = np.sum((cat0['WEAVE_priority1']) & (cat0['Total_flux']>=fluxcut) & (cat0['Maj']>size_large))
N13 = np.sum((cat13['WEAVE_priority1']) & (cat13['Total_flux']>=fluxcut) & (cat13['Maj']>size_large))
print ('0hr: {0:d}'.format(N0))
print ('13hr: {0:d}'.format(N13))

print()
print('WEAVE_priority 2')
N13 = np.sum((cat13['WEAVE_priority2']) & (cat13['Total_flux']>=fluxcut) & (cat13['Maj']>size_large))
print ('13hr: {0:d}'.format(N13))


print()
print('WEAVE_priority 3')
N13 = np.sum((cat13['WEAVE_priority3']) & (cat13['Total_flux']>=fluxcut) & (cat13['Maj']>size_large))
print ('13hr: {0:d}'.format(N13))

print()
print('HETDEX')
N13 = np.sum((cat13['HETDEX']) & (cat13['Total_flux']>=fluxcut) & (cat13['Maj']>size_large))
print ('13hr: {0:d}'.format(N13))


#print()
#print('WEAVE_priority 1')
#N0 = np.sum((cat0['WEAVE_priority1']) & (cat0['ID_flag']==5))
#N13 = np.sum((cat13['WEAVE_priority1']) & (cat13['ID_flag']==5))
##N13d = np.sum((cat13d['WEAVE_priority1']) & (cat13d['ID_flag']==5))
#print ('0hr: {0:d}'.format(N0))
#print ('13hr: {0:d}'.format(N13))
##print ('13hr (d): {0:d}'.format(N13d))

#print()
#print('WEAVE_priority 2')
#N13 = np.sum((cat13['WEAVE_priority2']) & (cat13['ID_flag']==5))
##N13 = np.sum((cat13d['WEAVE_priority2']) & (cat13d['ID_flag']==5))
#print ('13hr: {0:d}'.format(N13))
##print ('13hr (d): {0:d}'.format(N13d))


#print()
#print('WEAVE_priority 3')
#N13 = np.sum((cat13['WEAVE_priority3']) & (cat13['ID_flag']==5))
#print ('13hr: {0:d}'.format(N13))

print('13h')
print ('Flowchart LGZ selection')
for wpri in ['WEAVE_priority1','WEAVE_priority2','WEAVE_priority3','HETDEX']:
    sel =  (cat13[wpri]) & (
    #sel =   (
            (cat13['FC_flag2'] == 5) | \
            (cat13['FC_flag2'] == 12) | \
            (cat13['FC_flag2'] == 15) | \
            (cat13['FC_flag2'] == 26) | \
            (cat13['FC_flag2'] == 20) )
    N13 = np.sum(sel)
    print (' {0} 13hr: {1:d}'.format(wpri, N13))


print ('Flowchart PF selection')
for wpri in ['WEAVE_priority1','WEAVE_priority2','WEAVE_priority3','HETDEX']:
    sel =  (cat13[wpri]) & (
    #sel =   (
            (cat13['FC_flag2'] == 6) | \
            (cat13['FC_flag2'] == 13) | \
            (cat13['FC_flag2'] == 16) | \
            (cat13['FC_flag2'] == 27) | \
            (cat13['FC_flag2'] == 21) )
    N13 = np.sum(sel)
    print (' {0} 13hr: {1:d}'.format(wpri, N13))

print ('Flowchart PF LGZ selection')
for wpri in ['WEAVE_priority1','WEAVE_priority2','WEAVE_priority3','HETDEX']:
    sel =  (cat13[wpri]) & (cat13['Prefilter'] == 1) & (
    #sel =   (
            (cat13['FC_flag2'] == 6) | \
            (cat13['FC_flag2'] == 13) | \
            (cat13['FC_flag2'] == 16) | \
            (cat13['FC_flag2'] == 27) | \
            (cat13['FC_flag2'] == 21) )
    N13 = np.sum(sel)
    print (' {0} 13hr: {1:d}'.format(wpri, N13))


print ('Flowchart TBD selection')
for wpri in ['WEAVE_priority1','WEAVE_priority2','WEAVE_priority3','HETDEX']:
    sel =  (cat13[wpri]) & (
    #sel =   (
            (cat13['FC_flag2'] == 4) | \
            (cat13['FC_flag2'] == 9) | \
            (cat13['FC_flag2'] == 14) | \
            (cat13['FC_flag2'] == 17) | \
            (cat13['FC_flag2'] == 28) )
    N13 = np.sum(sel)
    N13ml = np.sum(sel & (cat13['ML_flag']==1))
    N13lgz = np.sum(sel & (cat13['ML_flag']==0))
    print (' {0} 13hr: {1:d}'.format(wpri, N13))
    print (' {0} ml 13hr: {1:d}'.format(wpri, N13ml))
    print (' {0} lgz 13hr: {1:d}'.format(wpri, N13lgz))

print ('Flowchart TBD selection large, faint')
for wpri in ['WEAVE_priority1','WEAVE_priority2','WEAVE_priority3','HETDEX']:
    sel =  (cat13[wpri]) & (
    #sel =   (
            (cat13['FC_flag2'] == 4) )
    N13 = np.sum(sel)
    N13ml = np.sum(sel & (cat13['ML_flag']==1))
    N13lgz = np.sum(sel & (cat13['ML_flag']==0))
    print (' {0} 13hr: {1:d}'.format(wpri, N13))
    print (' {0} ml 13hr: {1:d}'.format(wpri, N13ml))
    print (' {0} lgz 13hr: {1:d}'.format(wpri, N13lgz))

sys.exit()
print()
print('0hr')
print ('Flowchart LGZ selection')
for wpri in ['WEAVE_priority1']:
    #sel =  (cat0[wpri]) & (
    sel =   (
        (cat0['FC_flag2'] == 5) | \
        (cat0['FC_flag2'] == 9) | \
        (cat0['FC_flag2'] == 13) | \
        (cat0['FC_flag2'] == 15) | \
        (cat0['FC_flag2'] == 16) | \
        (cat0['FC_flag2'] == 20)| \
        (cat0['FC_flag2'] == 26) )
    N13 = np.sum(sel)
    print (' {0} 0hr: {1:d}'.format(wpri, N13))


print ('Flowchart PF selection')
for wpri in ['WEAVE_priority1']:
    #sel =  (cat0[wpri]) & (
    sel =   (
        (cat0['FC_flag2'] == 6) | \
        (cat0['FC_flag2'] == 14) | \
        (cat0['FC_flag2'] == 21) | \
        (cat0['FC_flag2'] == 27))
    N13 = np.sum(sel)
    print (' {0} 0hr: {1:d}'.format(wpri, N13))


print ('Flowchart TBD selection')
for wpri in ['WEAVE_priority1']:
    #sel =  (cat0[wpri]) & (
    sel =   (
        (cat0['FC_flag2'] == 4) | \
        (cat0['FC_flag2'] == 10) | \
        (cat0['FC_flag2'] == 17) | \
        (cat0['FC_flag2'] == 28))
    N13 = np.sum(sel)
    N13ml = np.sum(sel & (cat0['ML_flag']==1))
    N13lgz = np.sum(sel & (cat0['ML_flag']==0))
    print (' {0} 0hr: {1:d}'.format(wpri, N13))
    print (' {0} ml 13hr: {1:d}'.format(wpri, N13ml))
    print (' {0} lgz 13hr: {1:d}'.format(wpri, N13lgz))

print ('Flowchart TBD selection - large, faint')
for wpri in ['WEAVE_priority1']:
    #sel =  (cat0[wpri]) & (
    sel =   (
        (cat0['FC_flag2'] == 4) )
    N13 = np.sum(sel)
    N13ml = np.sum(sel & (cat0['ML_flag']==1))
    N13lgz = np.sum(sel & (cat0['ML_flag']==0))
    print (' {0} 0hr: {1:d}'.format(wpri, N13))
    print (' {0} ml 13hr: {1:d}'.format(wpri, N13ml))
    print (' {0} lgz 13hr: {1:d}'.format(wpri, N13lgz))


''' ID_flags are
-99 not set
0 - no id
1 - LR
2 - large optical galaxy
3 - LGZ
4 - visual id / prefilter
5 - tbd
6 - deblend
7 - too zoomed in after prefilter
8 - Uncatalogued host after prefilter
'''

#t,i = np.unique(cat0['ID_flag'], return_counts=True)
#for tt,ii in zip(t,i): print(tt,ii)
