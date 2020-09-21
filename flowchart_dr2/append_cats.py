import numpy as np
from astropy.table import Table, join, Column, vstack


path = '/data2/wwilliams/projects/lofar_surveys/LoTSS-DR2-Feb2020/lr/'

for t in ['srl','gaus']:
    for h in ['13h']:#['0h','13h']:
        scat1 = Table.read(path+'LoTSS_DR2_v100.{t}_s{h}.lr.fits'.format(h=h,t=t))
        scat2 = Table.read(path+'LoTSS_DR2_v100.{t}_n{h}.lr.fits'.format(h=h,t=t))
        scat1r = Table.read(path+'LoTSS_DR2_v100.{t}_s{h}.lr-remaining.fits'.format(h=h,t=t))
        scat2r = Table.read(path+'LoTSS_DR2_v100.{t}_n{h}.lr-remaining.fits'.format(h=h,t=t))
        
        scat = vstack((scat1,scat2,scat1r,scat2r))
        
        scat.sort('Source_Name')
        
        print('LoTSS_DR2_v100.{t}_{h}.lr.fits: '.format(h=h,t=t), len(scat), 'sources')
        scat.write(path+'LoTSS_DR2_v100.{t}_{h}.lr-full.fits'.format(h=h,t=t), overwrite=True)
        
        
    for h in ['0h']:#['0h','13h']:
        scat1 = Table.read(path+'LoTSS_DR2_v100.{t}_s{h}.lr.fits'.format(h=h,t=t))
        #scat2 = Table.read(path+'LoTSS_DR2_v100.{t}_n{h}.lr.fits'.format(h=h,t=t))
        scat1r = Table.read(path+'LoTSS_DR2_v100.{t}_s{h}.lr-remaining.fits'.format(h=h,t=t))
        #scat2r = Table.read(path+'LoTSS_DR2_v100.{t}_n{h}.lr-remaining.fits'.format(h=h,t=t))
        
        scat = vstack((scat1,scat1r))
        
        scat.sort('Source_Name')
        
        print('LoTSS_DR2_v100.{t}_s{h}.lr.fits: '.format(h=h,t=t), len(scat), 'sources')
        scat.write(path+'LoTSS_DR2_v100.{t}_{h}.lr-full.fits'.format(h=h,t=t), overwrite=True)




        
