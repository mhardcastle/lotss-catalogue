
import sys
import os
import numpy as np
from astropy.table import Table, Column, vstack
import astropy.io.fits as pf


if len(sys.argv) == 1:
    print("Usage is : python get_intersecting_sources.py field_code ")
    print('E.g.: python get_intersecting_sources.py 0 ')
    sys.exit(1)

h = str(sys.argv[1])
if 'h' not in h:
    h+='h'
if h not in  ['0h','13h','n0h','n13h','s0h','s13h']:
    print('unknown field code (should be 0h or 13h)',h)
    sys.exit(1)

path = '/data2/wwilliams/projects/lofar_surveys/LoTSS-DR2-Feb2020/'
lofarcat_file_srt = path+'LoTSS_DR2_v100.srl_{h}.lr-full.presort.fits'.format(h=h)

cat = Table.read(lofarcat_file_srt)

nsplit = 10

ra0 = cat['RA'].min()
ra1 = cat['RA'].max()
dec0 = cat['DEC'].min()
dec1 = cat['DEC'].max()

pad = 1.02*cat['Maj'].max()/3600.

ralims = np.linspace(0.99*ra0,1.01*ra1,nsplit+1)
declims = np.linspace(0.99*dec0,1.01*dec1,nsplit+1)

rac = (ralims[:-1] + ralims[1:])/2.
decc = (declims[:-1] + declims[1:])/2.
dra = np.mean(ralims[1:] - ralims[:-1])/2.
ddec = np.mean(declims[1:] - declims[:-1])/2.

if 'block' not in cat.colnames:
    cat.add_column(Column(data=-1*np.ones(len(cat)),name='block'))
#cat.add_column(Column(data=np.zeros(len(cat)),name='blockpad'))

do_calc = True
if do_calc:
    i =  0
    for ra in rac:
        for dec in decc:
            selpad = (cat['RA'] > (ra-dra-pad)) & (cat['RA'] <= (ra+dra+pad)) & (cat['DEC'] > (dec-ddec-pad)) & (cat['DEC'] <= (dec+ddec+pad))
            sel = (cat['RA'] > (ra-dra)) & (cat['RA'] <= (ra+dra)) & (cat['DEC'] > (dec-ddec)) & (cat['DEC'] <= (dec+ddec))
            
            cat['block'][sel] = i
            #cat['blockpad'][selpad] = i
            
            print(i, 'RA:',ra-dra-pad, ra+dra+pad, 'DEC:',dec-ddec-pad, dec+ddec+pad, len(cat[selpad]),len(cat[sel]))
            
            lofarcat_file_sub = path+'LoTSS_DR2_v100.srl_{h}_subcat{i:03d}.lr-full.presort.fits'.format(h=h,i=i)
            
            if not os.path.isfile(lofarcat_file_sub):
                
                if len(cat[selpad]) > 0:
                    cat[selpad].write(lofarcat_file_sub, overwrite=True)
                
                
                    os.system('python flow_python/get_intersecting_sources_multi.py {h}_subcat{i:03d}'.format(h=h,i=i))
            
            i+=1
            
        
      
# collate the blocks - we keep only the sources in the block, not the padded block
i =  0
t = None
for ra in rac:
    for dec in decc:
        
        sel = (cat['RA'] > (ra-dra)) & (cat['RA'] <= (ra+dra)) & (cat['DEC'] > (dec-ddec)) & (cat['DEC'] <= (dec+ddec))
        cat['block'][sel] = i
        
        lofarcat_file_sub = path+'LoTSS_DR2_v100.srl_{h}_subcat{i:03d}.lr-full.presort.fits'.format(h=h,i=i)
        print(i)
        
        if not os.path.isfile(lofarcat_file_sub):
            i+=1
            continue
        with pf.open(lofarcat_file_sub) as hd:
            ti = Table(hd[1].data)
        ti = ti[ti['block']==i]   # keep only block sources
        if len(ti)>0:
            if t is not None:
                t = vstack((t,ti))
            else:
                t = ti
            
        
        i += 1
cat.sort('Source_Name')
t.sort('Source_Name')

if 'Enclosed' not in  cat.colnames:
    cat.add_column(t['Enclosed'])
    cat.add_column(t['Encloses'])
    cat.add_column(t['Intersects'])
else:
    cat['Enclosed'] = t['Enclosed']
    cat['Encloses'] = t['Encloses']
    cat['Intersects'] = t['Intersects']

cat.write(lofarcat_file_srt, overwrite=True)

# get rid of the duplicates in the padded areas

'subcat'


if 0:
    i =  0
    for ra in rac:
        for dec in decc:
            
            selpad = (cat['RA'] > (ra-dra-pad)) & (cat['RA'] <= (ra+dra+pad)) & (cat['DEC'] > (dec-ddec-pad)) & (cat['DEC'] <= (dec+ddec+pad))
            sel = (cat['RA'] > (ra-dra)) & (cat['RA'] <= (ra+dra)) & (cat['DEC'] > (dec-ddec)) & (cat['DEC'] <= (dec+ddec))
