from astropy.table import Table
import numpy as np
from tqdm import tqdm
import healpy as hp
import matplotlib.pyplot as plt
import os
import pickle
import glob
from scipy.interpolate import interp1d

NSIDE=16 # pixels to work with

g=glob.glob('sourcecount-*.npy')

median_bins=np.logspace(-3,1.0,100)
median_bc=(median_bins[:-1]+median_bins[1:])/2
interps=[]
for f in tqdm(g):
    with open(f,'rb') as infile:
        bc,bins,dhist,hist=pickle.load(infile)

    bw=bins[1:]-bins[:-1]
    bw=bw[:len(dhist)]
    dhist/=bw
    i=interp1d(bc,dhist,bounds_error=False)
    interps.append(i(median_bc))

interps=np.array(interps)
print(interps.shape)

median=np.nanmedian(interps,axis=0)

with open('median.npy','wb') as outfile:
    pickle.dump([median_bc,median],outfile)

plt.plot(median_bc,median)
plt.xscale('log')
plt.yscale('log')
plt.show()
