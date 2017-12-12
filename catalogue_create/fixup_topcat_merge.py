from astropy.table import Table

# run the following

# /soft/topcat/stilts tskymatch2 ra1=ID_ra dec1=ID_dec ra2=ra dec2=dec error=10 out=merge.fits find=best1 join=all1 in1=LOFAR_HBA_T1_DR1_merge_ID_v0.6.fits in2=/data/lofar/mjh/hetdex_ps1_allwise_photoz_v0.2.fits

# and then run this to clean up afterwards.

t=Table.read('merge.fits')

t['RA_1'].name='RA'
t['DEC_1'].name='DEC'

t['ID_ra']=t['ra_2']
t['ID_dec']=t['dec_2']

t.remove_columns(['ra_2','dec_2','GroupID','GroupSize','Separation','raMean','decMean','class','id','z_good'])
t.sort('RA')
t.write('LOFAR_HBA_T1_DR1_merge_ID_optical_v0.3.fits',overwrite=True)
