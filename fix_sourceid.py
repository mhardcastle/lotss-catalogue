from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys

t=Table.read(sys.argv[1])
ilt=[]
sc=SkyCoord(t['RA'],t['DEC'],frame='icrs')
strings=sc.to_string(style='hmsdms',sep='',precision=2)
for s in strings:
    ilt.append(str('ILTJ'+s).replace(' ','')[:-1])
t['Source_Name']=ilt
if 'Maj' not in t.columns:
    # insert dummy
    t['Maj']=60
    t['Min']=60
    t['PA']=0

t.write(sys.argv[2])
