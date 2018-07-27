from astropy.table import Table

t=Table.read('LOFAR_HBA_T1_DR1_merge_ID_optical_v1.2.fits')
psf=t['RA']<184.445
psf&=t['DEC']>55.000
psf&=t['DEC']<55.2245

t['PanSTARRS_missing']=psf

t.write('LOFAR_HBA_T1_DR1_merge_ID_optical_f_v1.2.fits',overwrite=True)
