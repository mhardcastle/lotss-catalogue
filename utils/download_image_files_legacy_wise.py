from __future__ import print_function
from download_image_files import LofarMaps,get_legacy,get_first,get_wise
from find_wise import WISE
from astropy.table import Table
import sys
import os

if __name__=='__main__':
    t=Table.read(sys.argv[1])
    #if 'ra' in t.colnames:
    #    t['ra'].name='RA'
    #if 'dec' in t.colnames:
    #    t['dec'].name='DEC'
    outfilename=sys.argv[1].replace('.fits','-list.txt')
    if not(os.path.isfile(outfilename)):
        startpoint=0
        outfile=open(outfilename,'w')
    else:
        startpoint=len(open(outfilename).readlines())
        outfile=open(outfilename,'a')

    lm=LofarMaps(stay_in_imagedir=True)
    w=WISE()
    # now work in downloads dir
    os.chdir('downloads')

    for r in t[startpoint:]:
        print(r['Source_Name'])
        print(r['RA'],r['DEC'])
        lofarname=lm.find(r['RA'],r['DEC'])
        legacyname=get_legacy(r['RA'],r['DEC'],bands='zrg')
        wisename=w.find_pos(r['RA'],r['DEC'])
        if wisename is None:
            wisename=get_wise(r['RA'],r['DEC'],1)
        print(r['Source_Name'],lofarname,legacyname,wisename,file=outfile)
        outfile.flush()
    outfile.close()

