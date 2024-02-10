# Apply the initial filters

from astropy.table import Table

# Initial size > 20 arcsec
# Flux > 5 mJy

t=Table.read('/beegfs/lofar/mjh/rgz/combined-release-v1.0.fits')

filt=t['Total_flux']>5
filt&=t['LAS']>20
filt&=t['Resolved']

t[filt].write('sizeflux-in.fits',overwrite=True)
