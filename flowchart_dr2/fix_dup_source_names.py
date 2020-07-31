import string
import sys
import os
import numpy as np
from astropy.table import Table, Column, vstack




path = '/data2/wwilliams/projects/lofar_surveys/LoTSS-DR2-Feb2020/'

scat = Table.read (path+'lr/LoTSS_DR2_v100.{t}_{h}.lr-full.fits'.format(h='13h',t='srl'))
gcat = Table.read (path+'lr/LoTSS_DR2_v100.{t}_{h}.lr-full.fits'.format(h='13h',t='gaus'))

h='13h'
lofarcat_file_srt = path+'LoTSS_DR2_v100.srl_{h}.lr-full.presort.fits'.format(h=h)
scat = Table.read(lofarcat_file_srt)

# manually down for now ... but could be determined from np.unique
dup_names = ['ILTJ104642.23+371846.1',
'ILTJ113308.22+605332.4',
'ILTJ123019.45+564340.7',
'ILTJ125313.61+384929.8',
'ILTJ130927.32+500801.1',
'ILTJ163151.37+551815.7']

for dup in dup_names:
    si = np.where(scat['Source_Name'] == dup)[0]
    gi = np.where(gcat['Source_Name'] == dup)[0]
    
    if np.all(scat['S_Code'][si] =='S'):
        i = 0
        for sii in si:
        
            print(scat['Source_Name'][sii])
            print('clean')
            mgi = (gcat['RA'][gi] == scat['RA'][sii]) & (gcat['DEC'][gi] == scat['DEC'][sii]) & (gcat['Total_flux'][gi] == scat['Total_flux'][sii])
            scat['Source_Name'][sii] = scat['Source_Name'][sii]+string.ascii_letters[i]
            gcat['Source_Name'][gi[mgi]] = scat['Source_Name'][sii]
            
            i += 1
    else:
        
        # multiple gaus sources are painful
        # match them by finding the combination of gaus components that match their total flux
        # note no check is done that the final decomposition is mutually exclusive and spans the full gi set
        print(scat['Source_Name', 'S_Code'][si])
        print('messy gaus')
        
        allmgi = []
        i =0 
        for sii in si:
            sname =scat['Source_Name'][sii]
            flux = scat['Total_flux'][sii]
            print ('Source is',sname)
            print ('Source is',flux)
            gfluxes = gcat['Total_flux'][gi]
            
            
            
            Nm = len(gfluxes)
            
            # get all the combinations of gaus components - sinds is a list of unique arrays of binary indicies to gi for all the combnations
            s1 = [ bin(ll) for ll in range(2**Nm)]
            s1 = [ ss[2:] for ss in s1 ] 
            s1 = [ '0'*(Nm-len(ss))+ss  for ss in s1 ]
            sinds  = [ np.array([bool(int(ss)) for ss in si1 ]) for si1 in s1]
            
            # get sum of all combinations of components
            gfsum = np.array([ np.sum(gfluxes[sind]) for sind in sinds])
            
            # doesn't match exactly because of floating point representation errors
            sel = sinds[np.where(np.abs(gfsum - flux) < 1e-6)[0][0]]
            
            mgi = gi[sel]
            
            #print(mgi)
            
            scat['Source_Name'][sii] = scat['Source_Name'][sii]+string.ascii_letters[i]
            gcat['Source_Name'][mgi] = scat['Source_Name'][sii]
            
            # appedn mgi to  allmgi - will enable a check that all the mgi's were not repeated and al' the gi's were included
            
            i += 1
    
    print(scat[si])
    print(gcat[gi])


scat.write (path+'LoTSS_DR2_v100.{t}_{h}.lr-full.fits'.format(h='13h',t='srl'),overwrite=True)
gcat.write (path+'LoTSS_DR2_v100.{t}_{h}.lr-full.fits'.format(h='13h',t='gaus'),overwrite=True)
