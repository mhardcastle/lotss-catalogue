from astropy.table import Table,Column
import numpy as np
import sys
import os

# run the following

# /soft/topcat/stilts tskymatch2 ra1=ID_ra dec1=ID_dec ra2=ra dec2=dec error=5 out=merge.fits find=best1 join=all1 in1=../blend/merge_out.fits in2=/data/lofar/mjh/hetdex_ps1_allwise_photoz_v0.5.fits

# and then run this to clean up afterwards.

print 'Reading table'
t=Table.read('merge.fits')

print 'Rename columns'
t['RA_1'].name='RA'
t['DEC_1'].name='DEC'
t['ID_ra']=t['ra_2']
t['ID_dec']=t['dec_2']

print 'Remove unnecessary columns'
t.remove_columns(['ra_2','dec_2','Separation','raMean','decMean','class','id','z_good','Number_Masked','Number_Pointings'])

print 'Update ID names'
for i,r in enumerate(t):
    if ~np.isnan(r['ID_ra']) and r['ID_flag']!=2:
        psoname=r['objName'].rstrip()
        if psoname!="" and psoname!="N/A":
            t[i]['ID_name']=psoname
        else:
            t[i]['ID_name']='AllWISE '+r['AllWISE']

print 'Remove whitespace padding:',
sys.stdout.flush()
stringcols=[n for (n,ty) in t.dtype.descr if 'S' in ty]
for c in stringcols:
    print c,
    sys.stdout.flush()
    t[c]=[s.rstrip() for s in t[c]]
print

print 'Remove -99s:',
sys.stdout.flush()
dblcols=[n for (n,ty) in t.dtype.descr if 'f8' in ty]
for c in dblcols:
    print c,
    t[c]=np.where(t[c]==-99,np.nan,t[c])
print


print 'Sorting'
t.sort('RA')

print 'Writing to disk'
t.write('LOFAR_HBA_T1_DR1_merge_ID_optical_v1.1.fits',overwrite=True)
os.system('cp ../blend/merge_comp_out.fits LOFAR_HBA_T1_DR1_merge_ID_v1.1.comp.fits')
