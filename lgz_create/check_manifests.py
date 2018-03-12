#!/usr/bin/python

# Given a list file, check that the manifests exist and report if they don't.

import sys
import os

lines=open(sys.argv[1]).readlines()

for i,l in enumerate(lines):
    bits=l.split()
    if os.path.isfile(bits[0]+'-manifest.txt'):
        continue
    print i,bits[0]
