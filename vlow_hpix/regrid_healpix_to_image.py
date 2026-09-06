from astropy.io import fits
from reproject import reproject_from_healpix

hdu=fits.open('merged-weights-level-pixweight-4096.fits')
#hdu=fits.open('temperature_hp-8192.fits')
hdu[1].header['COORDSYS']='GALACTIC'

target_header = fits.Header.fromstring("""
NAXIS   =                    2
NAXIS1  =                  6000
NAXIS2  =                  2000
CTYPE1  = 'GLON-TAN'
CRPIX1  =              3000.5
CRVAL1  =                52.5
CDELT1  =               -0.01
CUNIT1  = 'deg     '
CTYPE2  = 'GLAT-TAN'
CRPIX2  =               1000.5
CRVAL2  =                  0.0
CDELT2  =                0.01
CUNIT2  = 'deg     '
COORDSYS= 'GALACTIC'
""", sep='\n')

array, footprint = reproject_from_healpix(hdu[1],target_header)
hdu_out=fits.PrimaryHDU(array,target_header)
#hdu_out.writeto('temp-8192-reproj.fits',overwrite=True)
hdu_out.writeto('feathered-4096-reproj.fits',overwrite=True)
