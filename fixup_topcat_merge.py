from astropy.table import Table

# run a topcat or stilts merge with max sep 10 arcsec on ID_ra and ra, ID_dec
# and dec. Save the result as merge.fits. Then

t=Table.read('merge.fits')

t['RA_1'].name='RA'
t['DEC_1'].name='DEC'

t['ID_ra']=t['ra_2']
t['ID_dec']=t['dec_2']

t.remove_columns(['ra_2','dec_2','GroupID','GroupSize','Separation','raMean','decMean','class','id'])
t.sort('RA')
t.write('LOFAR_HBA_T1_DR1_merge_ID_optical_v0.2.fits',overwrite=True)
