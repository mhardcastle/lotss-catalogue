#!/usr/bin/python

'''
make_summary
produce a summary pdfs with all the postage stamps in the given classes
'''


import glob
import os
import sys

# this should be the directory where the sample fits files are
ddir = sys.argv[1]

#look for all the sample fits files
fitslist = glob.glob(ddir+'/*.fits')

for f in fitslist:
    fdir = f.replace('.fits','')
    # get all the images for this sample
    filelist = glob.glob(fdir+'/*.png')
    print filelist
    cmd = 'montage '+' '.join(filelist)+'  -tile 5x -geometry 256x256+1+1 '+fdir+'.pdf'
    os.system(cmd)

