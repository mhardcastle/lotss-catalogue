from __future__ import print_function
from astropy.table import Table
import glob
import os
import numpy as np

mgt=Table.read('/beegfs/lofar/mjh/rgz/Spring/edited_gaussians.fits')

for field in ['Spring-40-45','Spring-60-65','Winter']:
    print(field)
    gt=Table.read('/beegfs/lofar/mjh/rgz/'+field+'/edited_gaussians.fits')
    g=glob.glob('/beegfs/lofar/mjh/rgz/'+field+'/blend/*.txt')
    for f in g:
        source=os.path.basename(f).replace('.txt','')
        print('        ',source)
        lines=open(f).readlines()
        incomps=False
        for i,l in enumerate(lines):
            if l.rstrip()=="":
                incomps=False
            if incomps:
                bits=l.split()
                index=int(bits[0])
                sname=gt[index]['Source_Name']
                findg=(mgt['Source_Name']==sname)
                assert(np.any(findg))
                newindex=np.argmax(findg)
                lines[i]=str(newindex)+' '+bits[1]+'\n'
            if l.startswith('## Components'):
                incomps=True
        with open('/beegfs/lofar/mjh/rgz/Spring/blend/'+source+'.txt','w') as newfile:
            newfile.write(''.join(lines))
            

            
    
