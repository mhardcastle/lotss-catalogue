import glob
import os
import emcee
import pickle
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from scipy.interpolate import interp1d

NSIDE=16

def loglf(cbins,hist,err,scale):
    mu=interp(cbins/10**scale)*10**scale
    likelihood=(hist-mu)**2/err**2
    return -0.5*np.nansum(likelihood)
    
def lnprior(X):
    if X[0]<-1 or X[0]>1:
        return -np.inf
    return 1

def lnpost(X,x,y,err):
    return lnprior(X)+loglf(x,y,err,X[0])

with open('median.npy','rb') as infile:
    median_bc,median=pickle.load(infile)

interp=interp1d(median_bc,median,bounds_error=False)

g=glob.glob('sourcecount-*.npy')
hmap=np.ones(hp.nside2npix(NSIDE))*np.nan
vals=[]
for f in g:
    print(f)
    picklefile=f
    pixel=int(f.replace('sourcecount-','').replace('.npy',''))
    with open(picklefile,'rb') as infile:
        bc,bins,dhist,hist=pickle.load(infile)
    bs=bins[1:]-bins[:-1]
    bs=bs[:len(dhist)]   
    hist=hist[:len(dhist)]
    dhist/=bs
    err=dhist/np.sqrt(hist)

    nwalkers=10
    ndim=1
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnpost,
                                    args=(bc,dhist,err))

    pos=[[0]+0.05*np.random.normal(size=ndim)
         for i in range(nwalkers)]

    sampler.run_mcmc(pos, 1000)
    samples=sampler.chain[:, 200:, :].reshape((-1, ndim))

    samplest=samples.transpose()

    means=np.mean(samplest,axis=1)
    errors=np.percentile(samplest,(16,84),axis=1)-means

    for i in range(ndim):
        print(i,means[i],errors[0,i],errors[1,i])

    fnorm=means[0]
    vals.append(10**fnorm)
    hmap[pixel]=10**fnorm

hp.fitsfunc.write_map('mscnorm_hp.fits',hmap,overwrite=True)
    
plt.hist(vals,range=(0.8,1.2),bins=50)
print('Mean is',np.mean(vals))
print('Median is',np.median(vals))
print(np.percentile(vals,(16,84)))
plt.show()
