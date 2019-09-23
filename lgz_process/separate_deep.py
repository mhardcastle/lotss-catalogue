# Take the output of galzoo_export.py and split it into three sets of CSV files
import glob
import sys
from astropy.coordinates import SkyCoord,get_icrs_coordinates
import astropy.units as u

g=glob.glob('lofar-*.csv')
for f in g:
    lines=open(f).readlines()
    bits=lines[0].split(',')
    for i in range(len(bits)):
        if bits[i]=='source_name':
            break
    else:
        print 'Index not found!'
        sys.exit(1)
    index=i

    print f,'index is',index
    fields=['bootes','lockman','en1']
    rar=[(214.0,222.0),(156.0,168.0),(237.0,248.0)]
    files=[]
    for fd in fields:
        files.append(open(fd+'/'+f,'w'))
        files[-1].write(lines[0])
    for l in lines[1:]:
        bits=l.split(',')
        name=bits[index]
        s=name[4:]
        coord=s[0:2]+':'+s[2:4]+':'+s[4:9]+' '+s[9:12]+':'+s[12:14]+':'+s[14:]
        sc = SkyCoord(coord,unit=(u.hourangle,u.deg))
        ra=sc.ra.value
        dec=sc.dec.value
        print name,ra,dec
        for j,rnge in enumerate(rar):
            if ra>rnge[0] and ra<rnge[1]:
                fno=j
                break
        else:
            print 'Failed to find field'
            sys.exit(2)
        files[fno].write(l)
        
for f in files:
    f.close()
