from astropy.table import Table
import glob
import sys
import os

g=sorted(glob.glob('sources-v*.fits'))
infile=g[-1]
ts=Table.read(infile)

g=sorted(glob.glob('components-v*.fits'))
infile=g[-1]
tc=Table.read(infile)

sourcename=sys.argv[1]

r=ts[ts['Source_Name']==sourcename]
if len(r)>0:
    print 'In source list'
    print r
    r.write('temp-list.fits',overwrite=True)
    r=tc[tc['Parent_Source']==sourcename]
    print 'Children are'
    print r
    if len(sys.argv)>2:
        os.system('python /home/mjh/git/lotss-catalogue/lgz_create/make_overlays_EN1_inspect.py temp-list.fits 0')
        os.system('rm temp-list.fits')
else:
    r=ts[ts['Renamed_from']==sourcename]
    if len(r)>0:
        print 'Renamed'
        print r
    else:
        r=tc[tc['Component_Name']==sourcename]
        if len(r)>0:
            print 'In component list'
            print r
        else:
            r=tc[tc['Deblended_from']==sourcename]
            if len(r)>0:
                print 'Was deblended, children are'
                print r
                
