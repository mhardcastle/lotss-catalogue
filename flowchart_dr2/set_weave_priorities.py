import sys
import numpy as np
from astropy.table import Table, Column
import astropy.units as u    
import astropy.coordinates as ac
from shapely.geometry import Polygon, Point
import mocpy




def inHETDEX(ra,dec):
    inside = np.zeros(len(ra), dtype=bool)
    
    hvert = Table.read('/home/wwilliams/data2/projects/lofar_surveys/DR1/hetdex_vertices.fits')
    hra = hvert['HETDEX_RA'][0]
    hdec = hvert['HETDEX_DEC'][0]
    
    phetdex = Polygon(np.array((hra, hdec)).T)
    for i in range(len(ra)):
        inside[i] = phetdex.contains(Point(ra[i],dec[i]))
    
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
    if 'WEAVE_priority1' not in lofarcat.colnames:
        lofarcat.add_column(Column(data=np.zeros(len(lofarcat),dtype=bool), name='WEAVE_priority1'))
        lofarcat.add_column(Column(data=np.zeros(len(lofarcat),dtype=bool), name='WEAVE_priority1a'))
        lofarcat.add_column(Column(data=np.zeros(len(lofarcat),dtype=bool), name='WEAVE_priority2'))
        lofarcat.add_column(Column(data=np.zeros(len(lofarcat),dtype=bool), name='WEAVE_priority3'))
    else:
        lofarcat['WEAVE_priority1'][:] = False
        lofarcat['WEAVE_priority1a'][:] = False
        lofarcat['WEAVE_priority2'][:] = False
        lofarcat['WEAVE_priority3'][:] = False

    ind8hr = (lofarcat['RA']>= 110) & (lofarcat['RA']<= 135) & (lofarcat['DEC']>= 26) & (lofarcat['DEC']<= 41)
    #ind13hr45 = (lofarcat['RA']>= 130) & (lofarcat['RA']<= 280) & (lofarcat['DEC']>= 40) & (lofarcat['DEC']<= 45)
    ind13hr45 = (lofarcat['DEC']>= 39.8) & (lofarcat['DEC']<= 46.05)
    ind13hr60_1a = (lofarcat['RA']>= 120) & (lofarcat['RA']<= 240) & (lofarcat['DEC']>= 60) & (lofarcat['DEC']<= 65)
    ind13hr60 = (lofarcat['RA']>= 120) & (lofarcat['RA']<= 240) & (lofarcat['DEC']>= 59) & (lofarcat['DEC']<= 65.5)
    
    # 0hr start with all
    ind0hr = (lofarcat['RA']>= 0) & (lofarcat['RA']<= 360) & (lofarcat['DEC']>= -90) & (lofarcat['DEC']<= 90)
    
    ## remove dr1 area:   - temp take out - slow and small area
    #indhetdex = inHETDEX(lofarcat['RA'],lofarcat['DEC'])
    #lofarcat['WEAVE_priority0'][indhetdex] = True
    
    lofarcat['WEAVE_priority1a'][ind13hr60_1a] = True  # this is the initial bit done - later expanded to full priority 1 area
    
    lofarcat['WEAVE_priority1'][ind13hr60] = True
    
    lofarcat['WEAVE_priority2'][ind8hr] = True  
    lofarcat['WEAVE_priority3'][ind13hr45& (~ind8hr)] = True    # 45 dec strip overlaps with 8hr region, so exclude that
    
    #lofarcat['WEAVE_priority1'][indhetdex] = False
    #lofarcat['WEAVE_priority3'][indhetdex] = False
    
    print('{} sources out of {} with WEAVE_priority1 ({:.1f}%)'.format(np.sum(lofarcat['WEAVE_priority1']), Nlofarcat, 100.*np.sum(lofarcat['WEAVE_priority1'])/Nlofarcat))
    print('{} sources out of {} with WEAVE_priority2 ({:.1f}%)'.format(np.sum(lofarcat['WEAVE_priority2']), Nlofarcat, 100.*np.sum(lofarcat['WEAVE_priority2'])/Nlofarcat))
    print('{} sources out of {} with WEAVE_priority3 ({:.1f}%)'.format(np.sum(lofarcat['WEAVE_priority3']), Nlofarcat, 100.*np.sum(lofarcat['WEAVE_priority3'])/Nlofarcat))
    
    
    # priority area to complete IDs is same as 8hr weave priority
    if 'Complete_priority1' not in lofarcat.colnames:
        lofarcat.add_column(Column(data=np.zeros(len(lofarcat),dtype=bool), name='Complete_priority1'))
    else:
        lofarcat['Complete_priority1'][:] = False
        
    lofarcat['Complete_priority1'][ind8hr] = True
    
elif h =='0h':
    if 'WEAVE_priority1' not in lofarcat.colnames:
        lofarcat.add_column(Column(data=np.zeros(len(lofarcat),dtype=bool), name='WEAVE_priority1'))
    
    c = ac.SkyCoord(lofarcat['RA'],lofarcat['DEC'])
    legacy_moc_s = mocpy.MOC.from_fits(path+'/mocs/legacy-south.moc')
    dr2_moc = mocpy.MOC.from_fits(path+'/mocs/dr2-moc.moc') 
    
    ## do intersection in Aladin - gives area 1087 sq deg
    inlegacy = legacy_moc_s.contains(c.ra, c.dec)
    lofarcat['WEAVE_priority1'][inlegacy] = True

    
    print('{} sources out of {} with WEAVE_priority1 ({:.1f}%)'.format(np.sum(lofarcat['WEAVE_priority1']), Nlofarcat, 100.*np.sum(lofarcat['WEAVE_priority1'])/Nlofarcat))
    
    
    # priority area to complete IDs is same as 8hr weave priority
    if 'Complete_priority1' not in lofarcat.colnames:
        lofarcat.add_column(Column(data=np.zeros(len(lofarcat),dtype=bool), name='Complete_priority1'))
    else:
        lofarcat['Complete_priority1'][:] = False
        
    # priority area to complete IDs is pointing P006+31
    indP006p31 = lofarcat['Mosaic_ID']=='P006+31' 
    lofarcat['Complete_priority1'][indP006p31] = True
    
lofarcat.write(lofarcat_file_srt, overwrite=True, serialize_meta=True)

            
