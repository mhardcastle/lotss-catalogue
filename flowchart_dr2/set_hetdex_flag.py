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

def inHETDEXm(mid):
    inside = np.zeros(len(ra), dtype=bool)
    hetdex_fields = ['P10Hetdex', 'P11Hetdex12', 'P12Hetdex11', 'P14Hetdex04', 'P15Hetdex13', 'P164+55', 'P169+55', 'P16Hetdex13', 'P173+55', 'P178+55', 'P182+55', 'P187+55', 'P18Hetdex03', 'P191+55', 'P196+55', 'P19Hetdex17', 'P1Hetdex15', 'P200+55', 'P205+55', 'P206+50', 'P206+52', 'P209+55', 'P21', 'P210+47', 'P211+50', 'P213+47', 'P214+55', 'P217+47', 'P218+55', 'P219+50', 'P219+52', 'P221+47', 'P223+50', 'P223+52', 'P223+55', 'P225+47', 'P227+50', 'P227+53', 'P22Hetdex04', 'P23Hetdex20', 'P25Hetdex09', 'P26Hetdex03', 'P27Hetdex09', 'P29Hetdex19', 'P30Hetdex06', 'P33Hetdex08', 'P34Hetdex06', 'P35Hetdex10', 'P37Hetdex15', 'P38Hetdex07', 'P39Hetdex19', 'P3Hetdex16', 'P41Hetdex', 'P42Hetdex07', 'P4Hetdex16', 'P6', 'P7Hetdex11', 'P8Hetdex']
    hetdex_fields.append(['P210+52', 'P214+52', 'P215+50', 'P31Hetdex19'])
        
    for i in range(len(mid)):
        inside[i] = mid[i] in hetdex_fields
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
    indhetdex = inHETDEXm(lofarcat['Mosaic_ID'])
    lofarcat['HETDEX'][indhetdex] = True
    
    print('{} sources out of {} with in hetdex ({:.1f}%)'.format(np.sum(lofarcat['WEAVE_priority1']), Nlofarcat, 100.*np.sum(lofarcat['HETDEX'])/Nlofarcat))
    
     
    
lofarcat.write(lofarcat_file_srt, overwrite=True, serialize_meta=True)

            
