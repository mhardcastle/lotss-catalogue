import os
import sys
import glob

name=sys.argv[1]

files=glob.glob('*/'+name+'-manifest.txt')
if len(files)==0:
    raise RuntimeError('Files do not exist')

dir=files[0].split('/')[0]
print 'Working in dir',dir
os.chdir(dir)    

flist=glob.glob('*-list.txt')

lines=open(flist[0]).readlines()
for i,l in enumerate(lines):
    if name in l:
        break
else:
    raise RuntimeError('Name '+name+' not found!')

print 'Found',name,'as list entry',i

os.system('rm '+name+'*')
fitsname=flist[0].replace('-list.txt','.fits')
os.system('python $LGZPATH/lgz_create/make_overlays_legacy_cscale.py %s %i' % (fitsname,i))
