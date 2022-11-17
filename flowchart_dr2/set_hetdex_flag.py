import sys
import numpy as np
from astropy.table import Table, Column
import astropy.units as u    
import astropy.coordinates as ac
from shapely.geometry import Polygon, Point
import mocpy




def inHETDEX(ra,dec):
    inside = np.zeros(len(ra), dtype=bool)
    
    hvert = Table.read('/Users/w.williams/projects/lofar_surveys/DR1/hetdex_vertices.fits')
    hra = hvert['HETDEX_RA'][0]
    hdec = hvert['HETDEX_DEC'][0]
    
    phetdex = Polygon(np.array((hra, hdec)).T)
    #phetdex = phetdex.convex_hull
    
    
    #hmoc = mocpy.MOC('/home/wwilliams/data2/projects/lofar_surveys/DR1/DR1_moc.fits')
    
    for i in range(len(ra)):
        inside[i] = phetdex.contains(Point(ra[i],dec[i]))
        #inside[i] = hmoc.contains(ra[i],dec[i])
    
    return inside

h = str(sys.argv[1])
if 'h' not in h:
    h+='h'
if h not in  ['0h','13h']:
    print('unknown field code (should be 0h or 13h)',h)
    sys.exit(1)

path = '/Users/w.williams/projects/lofar_surveys/DR2/'

version ='v110'

lofarcat_file = path+'LoTSS_DR2_{version}.srl_{h}.lr-full.presort.hdf5'.format(h=h,version=version)
lofarcat_file_srt = lofarcat_file
lofarcat = Table.read(lofarcat_file)


Nlofarcat = len(lofarcat)

if h =='13h':
    if 'HETDEX' not in lofarcat.colnames:
        lofarcat.add_column(Column(data=np.zeros(len(lofarcat),dtype=bool), name='HETDEX'))
    else:
        lofarcat['HETDEX'][:] = False

    
    ## remove dr1 area:   - temp take out - slow and small area
    indhetdex = inHETDEX(lofarcat['RA'],lofarcat['DEC'])
    lofarcat['HETDEX'][indhetdex] = True
    
    print('{} sources out of {} with in hetdex ({:.1f}%)'.format(np.sum(lofarcat['WEAVE_priority1']), Nlofarcat, 100.*np.sum(lofarcat['HETDEX'])/Nlofarcat))
    
     
    
lofarcat.write(lofarcat_file_srt, overwrite=True, serialize_meta=True)

            
