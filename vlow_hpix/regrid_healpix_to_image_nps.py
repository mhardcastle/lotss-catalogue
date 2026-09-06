from astropy.io import fits
from reproject import reproject_from_healpix

for i,image in enumerate(['merged-weights-level-pixweight-4096.fits','haslam408_dsds_Remazeilles2014.fits','ovro_lwa_sky_map_73.152MHz.fits']):
    
    hdu=fits.open(image)
    hdu[1].header['COORDSYS']='GALACTIC'

    target_header = fits.Header.fromstring("""
NAXIS   =                    2
NAXIS1  =                  1200
NAXIS2  =                  3000
CTYPE1  = 'GLON-TAN'
CRPIX1  =               600.5
CRVAL1  =                32.0
CDELT1  =               -0.01
CUNIT1  = 'deg     '
CTYPE2  = 'GLAT-TAN'
CRPIX2  =              2000.5
CRVAL2  =                34.2
CDELT2  =                0.01
CUNIT2  = 'deg     '
COORDSYS= 'GALACTIC'
""", sep='\n')

    array, footprint = reproject_from_healpix(hdu[1],target_header)
    if i==1:
        array*=(408/144)**2.7
    elif i==2:
        array*=(73/144)**2.7
    hdu_out=fits.PrimaryHDU(array,target_header)
    hdu_out.writeto(f'nps-{i}.fits',overwrite=True)
