#!/usr/bin/python

# Attempt to cross-match tables making use of the errors in a correct
# way.

from astropy.io import fits
#from astropy.io.votable import parse
import numpy as np
from numpy import rec
import matplotlib.pyplot as plt
from matplotlib import rc
from tqdm import tqdm
from astropy.table import Table

fontsize=16
rc('font',**{'family':'serif','serif':['Times'],'size':fontsize})
rc('text', usetex=True)

def getprop(table,field,lookuparray):
    n=len(lookuparray)
    temp=np.zeros(n)
    for i in range(0,n):
        if (lookuparray[i]>0):
            temp[i]=table[field][lookuparray[i]-1]
    return temp[lookuparray>0]

plt.figure(figsize=(6,6))

t_6c=Table.read('6C-DR3.fits')
n6c=len(t_6c)

snr_6c=t_6c['Peak_flux']/(0.13/5)
sigma_6c=np.sqrt(1.0+(74/snr_6c)**2.0)

print('Ingested',n6c,'6C sources')
print('Mean 6C positional error is',sigma_6c.mean(),'arcsec')

#plt.hist(sigma_6c,20)
#plt.show()

t_dr3=Table.read('bright-cut-NN1-6Cmatch.fits')
print(np.mean(t_dr3['Total_flux']))
t_dr3['Total_flux']/=1000
t_dr3['E_Total_flux']/=1000
t_dr3['Peak_flux']/=1000
t_dr3['E_Peak_flux']/=1000
print(np.mean(t_dr3['Total_flux']))
print(np.sum(t_dr3['Total_flux']>0.4))

ndr3=len(t_dr3)
snr_dr3=t_dr3['Peak_flux']/t_dr3['E_Peak_flux']
sigma_dr3=np.sqrt(0.5*(t_dr3['E_RA'][:]**2.0 + t_dr3['E_DEC'][:]**2.0))

print('Ingested',ndr3,'DR3 sources')
print('Mean DR3 positional error is',sigma_dr3.mean(),'arcsec')

plt.scatter(sigma_dr3,snr_dr3)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Positional error (arcsec)')
plt.ylabel('SNR')
plt.xlim((0.01,100))
r=np.logspace(-2,2)
condon=60.0/r/np.sqrt(8)
plt.plot(r,condon,color='red')
plt.savefig('snr.pdf')
plt.clf()

print('Filtering bad positions and NaN 151-MHz fluxes...')
filter=(sigma_dr3<20) & ~np.isnan(t_dr3['Peak_flux'])
t_dr3=t_dr3[filter]
ndr3=len(t_dr3)
sigma_dr3=sigma_dr3[filter]

print(ndr3,'DR3 sources after filtering')
print('Mean DR3 positional error after filtering is',sigma_dr3.mean(),'arcsec')

# Max likelihood cross-match

delta_ra=0.882 #7.50
delta_dec=0.279 # -0.668445223668 

ra=t_dr3['RA']
dec=t_dr3['DEC']
ra-=delta_ra/3600.0
dec-=delta_dec/3600.0

mmss=np.zeros(ndr3,dtype=int)
m6c=np.zeros(n6c,dtype=int)
l6c=np.zeros(ndr3)
off=np.zeros(ndr3)
for (i,r) in tqdm(list(enumerate(t_dr3))):
    ra=r['RA']
    dec=r['DEC']
    offset=np.sqrt((np.cos(dec*np.pi/180))*(t_6c['RA'][:]-ra)**2.0+(t_6c['DEC'][:]-dec)**2.0)*3600.0
    sigma2=sigma_dr3[i]**2.0+sigma_6c**2.0
    fsigma2=(0.1*r['Total_flux'])**2.0+(0.1*t_6c['Total_flux'])**2.0
    likelihood=np.log(offset/sigma2)-offset**2.0/(2.0*sigma2)-(r['Total_flux']-t_6c['Total_flux'])**2.0/(2.0*fsigma2)
    maxind=np.argmax(likelihood)
    if (likelihood[maxind]>-100.0):
        # check whether we are matching to a source that already has a match
        if (m6c[maxind]>0):
            j=m6c[maxind]-1
            if l6c[j]>likelihood[maxind]:
                continue
            else:
                #supersede that match
                mmss[j]=0
        mmss[i]=maxind+1
        m6c[maxind]=i+1
        l6c[i]=likelihood[maxind]
        off[i]=offset[maxind]

t_matched=t_dr3[mmss>0]
indices=mmss[mmss>0]-1
t6c_matched=t_6c[indices]
for c in t_6c.colnames:
    t_matched['6C_'+c]=t6c_matched[c]

t_matched['Ratio']=t_matched['Total_flux']/t_matched['6C_Total_flux']
t_matched.write('DR3-6C-MLM.fits')
        
#plt.hist(l6c,20,range=(-100,0))
#plt.scatter(-l6c,off)
#plt.show()

#plt.hist(off,20,range=(0,40))
#plt.show()

print(np.count_nonzero(mmss),'DR3 sources matched')
print(np.count_nonzero(m6c),'6C sources matched')

# Look at offsets

matched=t_dr3[mmss>0]
matchedmm=mmss[mmss>0]
dra=3600.0*(matched['RA']-getprop(t_6c,'RA',matchedmm))
ddec=3600.0*(matched['DEC']-getprop(t_6c,'DEC',matchedmm))

print('Mean delta RA is',dra.mean(),'and rms',dra.std())
print('Mean delta Dec is',ddec.mean(),'and rms',ddec.std())

print('Errors on the mean are',dra.std()/np.sqrt(len(matched)),ddec.std()/np.sqrt(len(matched)))

bright_dr3_filter=(t_dr3['Total_flux']>0.4)
bright_6c_filter=(t_6c['Total_flux']>0.4)
bmss=t_dr3[bright_dr3_filter]
b6c=t_6c[bright_6c_filter]
mbmss=mmss[bright_dr3_filter]
mb6c=m6c[bright_6c_filter]

print('Total number of DR3 sources with fluxes above 0.4 Jy is',len(bmss))
print('Total number of matched DR3 sources with fluxes above 0.4 Jy is',np.count_nonzero(mbmss))
print('Total number of 6C sources with fluxes above 0.4 Jy is',len(b6c))
print('Total number of matched 6C sources with fluxes above 0.4 Jy is',np.count_nonzero(mb6c))

matched=bmss[mbmss>0]
matchedmm=mbmss[mbmss>0]
dra=3600.0*(matched['RA']-getprop(t_6c,'RA',matchedmm))
ddec=3600.0*(matched['DEC']-getprop(t_6c,'DEC',matchedmm))

print('Mean delta RA is',dra.mean(),'and rms',dra.std())
print('Mean delta Dec is',ddec.mean(),'and rms',ddec.std())

print('Errors on the mean are',dra.std()/np.sqrt(len(matched)),ddec.std()/np.sqrt(len(matched)))

# plot positions of all matched sources

def plotblob(table,filt,ra,dec,flux,colour,label):
    rap=table[ra][filt]
    decp=table[dec][filt]
    #size=7*(2.4+np.log(table[flux][filt]))
    size=20
    plt.scatter(rap,decp,s=size,color=colour,label=label,alpha=0.6)
'''
#plt.ylim((63.5,76.5))
plotblob(t_dr3,mmss>0,'RA','DEC','Peak_flux','gray','Matched')
plotblob(t_dr3,(mmss==0),'RA','DEC','Peak_flux','blue','DR3 unmatched')
plotblob(t_6c,(m6c==0),'RA','DEC','Peak_flux','red','6C unmatched')

plt.xlabel('Right Ascension (J2000)')
plt.ylabel('Declination (J2000)')
plt.legend(loc=0,prop={'size':fontsize})

plt.savefig('sky-all.pdf')
plt.clf()
'''
# plot only the bright ones

#plt.ylim((63.5,77))
plotblob(bmss,mbmss>0,'RA','DEC','Peak_flux','green','Matched')
plotblob(bmss,(mbmss==0),'RA','DEC','Peak_flux','blue','DR3 unmatched')
plotblob(b6c,(mb6c==0),'RA','DEC','Peak_flux','red','6C unmatched')

plt.xlabel('Right Ascension (J2000)')
plt.ylabel('Declination (J2000)')
plt.legend(loc=0,prop={'size':fontsize})

plt.savefig('sky-bright.pdf')
plt.clf()
'''
# Print out the weird ones

for i in range(0,len(b6c)):
    if (mb6c[i]==0): print('6C',i,b6c[i]['RA'],b6c[i]['DEC'],b6c[i]['Peak_flux'])
print('---')
#for i in range(0,len(bmss)):
#    if (mbmss[i]==0): print('DR3',i,bmss[i]['Source_Name'],bmss[i]['RA'],bmss[i]['DEC'],bmss[i]['Peak_flux'])

# Total, bright sources only
'''
flpx=bmss['Total_flux'][(mbmss>0)]
flpxe=bmss['E_Total_flux'][(mbmss>0)]
flpy=getprop(t_6c,'Total_flux',mbmss)
flpye=0.1*flpy
plt.errorbar(flpx,flpy,xerr=flpxe,yerr=flpye,linestyle='none',alpha=0.1)
plt.plot([0.1,100],[0.1,100],color='red',linewidth=2.5)
plt.xlabel('DR3 total flux density (Jy)')
plt.ylabel('6C total flux density (Jy)')
plt.xscale('log')
plt.yscale('log')
plt.axis('equal')
plt.axis([0.1,100,0.1,100])

plt.savefig('ff-total.pdf')
plt.clf()

weights=1.0/((flpxe/flpx)**2.0+(flpye/flpy)**2.0)
print('Ratio of 6C to DR3 total fluxes is',np.average(flpy/flpx,weights=weights))

flpx=bmss['Peak_flux'][(mbmss>0)]
flpxe=bmss['E_Peak_flux'][(mbmss>0)]
flpy=getprop(t_6c,'Peak_flux',mbmss)
flpye=0.1*getprop(t_6c,'Peak_flux',mbmss)
plt.errorbar(flpx,flpy,xerr=flpxe,yerr=flpye,linestyle='none',alpha=0.1)
plt.plot([0.1,100],[0.1,100],color='red',linewidth=2.5)
plt.xlabel('DR3 peak flux density (Jy)')
plt.ylabel('6C peak flux density (Jy)')
plt.xscale('log')
plt.yscale('log')
plt.axis('equal')
plt.axis([0.1,100,0.1,100])

plt.savefig('ff-peak.pdf')
plt.clf()

weights=1.0/((flpxe/flpx)**2.0+(flpye/flpy)**2.0)
print('Ratio of 6C to DR3 peak fluxes is',np.average(flpy/flpx,weights=weights))

# Plot histograms of matched and unmatched source fluxes

_,bins,_=plt.hist(np.log10(t_dr3['Total_flux']),40,color='grey',label='All',alpha=0.1)
flpx=t_dr3['Total_flux'][(mmss>0)]
plt.hist(np.log10(flpx),bins=bins,color='green',label='Matched')
flpx=t_dr3['Total_flux'][(mmss==0)]
flpxf=flpx[~np.isnan(flpx)]
plt.hist(np.log10(flpxf),bins=bins,color='blue',alpha=0.3,label='DR3 unmatched')
flpx=t_6c['Total_flux'][m6c==0]
plt.hist(np.log10(flpx),bins=bins,color='red',alpha=0.3,label='6C unmatched')
plt.xlabel('log10(151-MHz flux density/Jy)')
plt.legend(loc=0,prop={'size':fontsize})
plt.plot([-0.4,-0.4],[0,120],color='yellow',linewidth=1)
plt.savefig('dethist.pdf')
plt.clf()

