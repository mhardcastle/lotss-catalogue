from __future__ import print_function
from download_image_files import LofarMaps,get_legacy,get_first,get_wise
from astropy.table import Table
import sys
import os

if __name__=='__main__':
    t=Table.read(sys.argv[1])
    if 'Dec' in t.colnames:
        t['Dec'].name='DEC'
    outfilename=sys.argv[1].replace('.fits','-list.txt')
    if not(os.path.isfile(outfilename)):
        startpoint=0
        outfile=open(outfilename,'w')
    else:
        startpoint=len(open(outfilename).readlines())
        outfile=open(outfilename,'a')

    lm=LofarMaps(stay_in_imagedir=True)
    
    # now work in downloads dir
    os.chdir('downloads')

    for r in t[startpoint:]:
        print(r['Source_Name'])
        print(r['RA'],r['DEC'])
        lofarname=lm.find(r['RA'],r['DEC'])
        try:
            legacyname=get_legacy(r['RA'],r['DEC'],bands='zrg')
        except RuntimeError as e: # assume catching a 500 here
            if '500' in e.args[0]:
                print('Caught a 500, no image!')
                legacyname="None"
            else:
                raise
                
        #wisename=get_wise(r['RA'],r['DEC'],1)
        #firstname=get_first(r['RA'],r['DEC'])
        print(r['Source_Name'],lofarname,legacyname,file=outfile) 

    outfile.close()

