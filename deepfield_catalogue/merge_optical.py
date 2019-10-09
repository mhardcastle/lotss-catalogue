import glob
import os
from filter_cat import filter_table
from astropy.table import Table
from final_flag_clean_overlap import final_flag
import sys
import numpy as np

dir=os.getcwd()
field=os.path.basename(dir)
print 'field is',field

if field=='en1':
    mask='/beegfs/lofar/deepfields/ELAIS_N1_optical/radio_optical_overlap_masks/image_full_ampphase_di_m.NS_shift.blanked.scaled.rms_spmask.fits'
    optcat='/beegfs/lofar/deepfields/ELAIS_N1_optical/catalogues/correct_merging/add_uncat/EN1_MASTER_opt_spitzer_merged_cedit_apcorr_adduncat_lite.fits'
    src='/beegfs/lofar/deepfields/science_ready_catalogs/EN1_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits'
else:
    raise RuntimeError('Field not supported!')
    
g=sorted(glob.glob('sources-v*.fits'))

infile=g[-1]
print 'Processing',infile
outfile=infile.replace('sources','filtered')

filter_table(infile,mask,outfile)

flagname=outfile.replace('filtered','flagged')
print 'Applying final flags'
final_flag(field,outfile,flagname)

mergeout=flagname.replace('flagged','merged')

os.system('/soft/topcat/stilts tskymatch2 ra1=optRA dec1=optDec ra2=ALPHA_J2000 dec2=DELTA_J2000 error=1.5 join=all1 in1=%s in2=%s out=%s' % (flagname,optcat,mergeout))

mergeout2=mergeout.replace('merged','merged_src')

os.system('/soft/topcat/stilts tmatch2  join=all1 values1=NUMBER values2=ID matcher=exact in1=%s in2=%s out=%s' % (mergeout,src,mergeout2))

t=Table.read(mergeout2)
finalname=mergeout2.replace('merged_src','final')
print 'Remove unnecessary columns'
t.remove_columns(['RA_2','DEC_2','FLAG_OVERLAP_2','FLAG_CLEAN_2','id', 'ID_OPTICAL', 'ID_SPITZER'])
print 'Rename columns'
#t['Separation_1'].name='Separation'
for column in ['RA','DEC','FLAG_OVERLAP','flag_clean']:
    t[column+'_1'].name=column

print 'Remove whitespace padding:',
sys.stdout.flush()
stringcols=[n for (n,ty) in t.dtype.descr if 'S' in ty]
for c in stringcols:
    print c,
    sys.stdout.flush()
    t[c]=[s.rstrip() for s in t[c]]
print

print 'Remove -99s:',
dblcols=[n for (n,ty) in t.dtype.descr if ('f8' in ty or 'f4' in ty)]
for c in dblcols:
    print c,
    sys.stdout.flush()
    t[c]=np.where(t[c]==-99,np.nan,t[c])
print


print 'Sorting'
t.sort('RA')

print 'Writing to disk'
t.write(finalname,overwrite=True)

