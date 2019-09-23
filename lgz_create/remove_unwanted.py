import os
from astropy.table import Table

t_old=Table.read('lgz.fits')
t_new=Table.read('lgz-v2.fits')

os.system('mkdir unwanted')

for r in t_old:
    name=r['Source_Name']
    if name not in t_new['Source_Name']:
        print 'Removing',name
        os.system('mv '+name+'* unwanted')

