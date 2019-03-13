import glob
from PIL import Image
import numpy as np

outfile=open('manifest.csv','w')
outfile.write('subject_id,image_name_1,image_name_2,image_name_3,source_name,ra,dec,#size\n')
g=sorted(glob.glob('*-manifest.txt'))
for f in g:
    source=f.replace('-manifest.txt','')
    pngs=sorted(glob.glob(source+'*.png'))
    print source,
    sizes=[]
    for p in pngs:
        with Image.open(p) as img:
            print img.size,
            sizes.append(img.size)
    sizes=np.array(sizes)
    xvar=np.max(sizes[:,0])-np.min(sizes[:,0])
    yvar=np.max(sizes[:,1])-np.min(sizes[:,1])
    if xvar/np.mean(sizes[:,0]) > 0.05:
        problem=True
    elif yvar/np.mean(sizes[:,1]) > 0.05:
        problem=True
    else:
        problem=False
    if problem:
        print 'problem',
    print
    if not problem:
        line=open(f).readlines()[0]
        outfile.write(line)
