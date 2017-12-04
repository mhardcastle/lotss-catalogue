#!/usr/bin/python

'''
get_large_galaxies
produce fits catalogues of all the lofar sources corresponding to known nearby galaxies (based on their 2MASX names
'''

############################################################

clobber_2MAScat = True

path = '/local/wwilliams/projects/radio_imaging/lofar_surveys/LoTSS-DR1-July21-2017/'
lofarcat_file_srt = path+'LOFAR_HBA_T1_DR1_catalog_v0.9.srl.sorted.fits'


############################################################

from astropy.table import Table
import os

lofarcat_srt = Table.read(lofarcat_file_srt)

# Large nearby galaxies - 2MASX > 300"
filter_large = [['14031258+5420555','NGC5457'],
                ['12185761+4718133','M106', 'NGC4258'],
                ['13295269+4711429','M51a', 'NGC5194'],
                ['15155368+5619438','NGC5907'],
                ['11113096+5540268','M108', 'NGC3556']]


for tl in filter_large:
    t = tl[0]
    tn = tl[1]
    ind = lofarcat_srt['2MASX_name'] == t
    ll = lofarcat_srt[ind]
    catname = 'LOFAR_{0:s}.fits'.format(tn)
    if os.path.exists(catname):
        if clobber_2MAScat:
            os.remove(catname)
            ll.write(catname)
        else:
            print 'cannot write ',catname
            
    else:
        ll.write(catname)