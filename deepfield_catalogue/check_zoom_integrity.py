from __future__ import print_function
import glob

g=glob.glob('zoom/*.txt')

d={}

for f in g:
    z=[l.rstrip() for l in open(f).readlines()]
    for l in z:
        if l.startswith('ILTJ'):
            if l in d:
                print('Zoom file',f,'has duplicate source',l)
                print('Also in',d[l])
            else:
                d[l]=f
                
