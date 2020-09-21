from __future__ import print_function, division
from astropy.table import Table
import sys
import numpy as np
from surveys_db import SurveysDB
import matplotlib.pyplot as plt

p=[]
f=[]
b=[]
z=[]
a=[]

with SurveysDB(readonly=True) as sdb:

    filenames=sys.argv[1:]

    for filename in filenames:
        field=filename.split('/')[0]
        r=sdb.get_field(field)
        if r is not None:
            priority=r['weave_priority']
        else:
            priority=-1
        t=Table.read(filename)
        s=len(t)

        print('\n---------\n%s (%i)' % (filename,priority))
        print('Total number of sources is',s)

        print('Total fraction with opt ID is',np.sum(~np.isnan(t['optRA']))/s)
        print('Total fraction of associations',np.sum(t['Assoc']>0)/s)
        print('Total fraction with TZI is',np.sum(t['Zoom_prob']>0.5)/s)
        print('Total fraction with Blend is',np.sum(t['Blend_prob']>0.5)/s)
        print('Total fraction with Artefact is',np.sum(t['Art_prob']>0.5)/s)
        print('Total fraction with Imagemissing is',np.sum(t['Imagemissing_prob']>0.6)/s)
        print('Total fraction with Badclick>1 is',np.sum(t['Badclick']>1)/s)
        p.append(priority)
        f.append(np.sum(~np.isnan(t['optRA']))/s)
        z.append(np.sum(t['Zoom_prob']>0.5)/s)
        b.append(np.sum(t['Badclick']>1)/s)
        a.append(np.sum(t['Assoc']>0)/s)

plt.subplot(2,2,1)
plt.scatter(p[:-1],f[:-1])
print('Mean ID fraction',np.mean(f[:-1]))
plt.scatter(60,f[-1])
plt.ylabel('Fraction of majority IDs')

plt.subplot(2,2,2)
plt.scatter(p[:-1],z[:-1])
print('Mean TZI',np.mean(z[:-1]))
plt.scatter(60,z[-1])
plt.ylabel('Fraction of TZI')

plt.subplot(2,2,3)
plt.scatter(p[:-1],b[:-1])
print('Mean badclick',np.mean(b[:-1]))
plt.scatter(60,b[-1])
plt.xlabel('Priority (order)')
plt.ylabel('Fraction of badclick>1')

plt.subplot(2,2,4)
plt.scatter(p[:-1],a[:-1])
print('Mean assoc',np.mean(a[:-1]))
plt.scatter(60,a[-1])
plt.xlabel('Priority (order)')
plt.ylabel('Fraction of associations')

plt.show()
