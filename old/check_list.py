import os
import sys
lines=open(sys.argv[1]).readlines()

names=[l.split()[0] for l in lines]

for i,n in enumerate(names):
    if not os.path.isfile(n+'-manifest.txt'):
        print i,n
