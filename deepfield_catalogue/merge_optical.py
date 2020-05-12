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
try:
    skip_noid=(sys.argv[1]=='noid')
except:
    skip_noid=False

if field=='en1':
    mask='/beegfs/lofar/deepfields/ELAIS_N1_optical/radio_optical_overlap_masks/image_full_ampphase_di_m.NS_shift.int.facetRestored-scaled.pybdsm.rmsd_I_spmask.fits'
    optcat='/beegfs/lofar/deepfields/ELAIS_N1_optical/catalogues/correct_merging/add_uncat/EN1_MASTER_opt_spitzer_merged_cedit_apcorr_adduncat_lite.fits'
    src='/beegfs/lofar/deepfields/science_ready_catalogs/EN1_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_masses_public.fits'
    xid='/beegfs/lofar/deepfields/ELAIS_N1_FIR_prelim/XID+_lofar_ELAIS-N1_v0.5_20200113.fits'
    src_ze='/beegfs/lofar/deepfields/science_ready_catalogs/EN1_opt_spitzer_merged_vac_opt3as_irac4as_specz.v2.fits'
    new_lr_path = "/beegfs/lofar/deepfields/LR_update/en1/EN1_ML_RUN_fin_overlap_srl.fits"
    lr_vis_path = "/beegfs/lofar/deepfields/LR_update/en1/EN1_lr_in_old_not_new_pos.txt"
    lr_th = 0.056
    changeid_path = "/beegfs/lofar/deepfields/LR_update/en1/EN1_change_IDs.txt"
elif field=='bootes':
    mask='/beegfs/lofar/deepfields/Bootes_merged_optical/radio_optical_overlap_masks/image_full_ampphase_di_m.NS_shift.blanked.scaled.rms_spmask.fits'
    optcat='/beegfs/lofar/deepfields/Bootes_merged_optical/add_uncat/Bootes_MASTER_opt_spitzer_merged_adduncat_lite.fits'
    src='/beegfs/lofar/deepfields/science_ready_catalogs/Bootes_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_masses_public.fits'
    xid='/beegfs/lofar/deepfields/Bootes_FIR/XID+_lofar_Bootes_v0.5_20200209.fits'
    src_ze=None
    new_lr_path = "/beegfs/lofar/deepfields/LR_update/bootes/Bootes_ML_RUN_fin_overlap_srl.fits"
    lr_vis_path = "/beegfs/lofar/deepfields/LR_update/bootes/Bootes_lr_in_old_not_new_pos.txt"
    lr_th = 0.22
    changeid_path = None
elif field=='lockman':
    mask='/beegfs/lofar/deepfields/Lockman_edited_cats/radio_optical_overlap_masks/image_full_ampphase_di_m.NS_shift.blanked.scaled.rms_spmask.fits'
    optcat='/beegfs/lofar/deepfields/Lockman_edited_cats/optical/add_uncat/LH_MASTER_opt_spitzer_merged_cedit_apcorr_adduncat_lite.fits'
    src='/beegfs/lofar/deepfields/science_ready_catalogs/LH_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_masses_public.fits'
    xid='/beegfs/lofar/deepfields/Lockman_FIR/XID+_lofar_Lockman_v0.5_20200303.fits'
    src_ze='/beegfs/lofar/deepfields/science_ready_catalogs/LH_opt_spitzer_merged_vac_opt3as_irac4as_specz.v2.fits'
    new_lr_path = "/beegfs/lofar/deepfields/LR_update/lockman/LH_ML_RUN_fin_overlap_srl.fits"
    lr_vis_path = "/beegfs/lofar/deepfields/LR_update/lockman/LH_lr_in_old_not_new_pos.txt"
    lr_th = 0.055
    changeid_path = "/beegfs/lofar/deepfields/LR_update/lockman/Lockman_change_IDs.txt"
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

#if field=='bootes':
#    print 'Fixing ID column!'
#    t=Table.read(mergeout)
#    t['NUMBER']=t['ID'].astype('int')
#    t['NUMBER'].mask=np.isnan(t['ID'])
#    t.write(mergeout,overwrite=True)

if src_ze is not None:
    # Join the src and src_ze columns
    src_t = Table.read(src)
    src_ze_t = Table.read(src_ze)

    new_z_cols = ["Z_SPEC", "Z_SOURCE", "Z_QUAL"]
    for col in new_z_cols:
        src_t[col] = src_ze_t[col]
    del src_ze_t
    # Write this to a tmp file
    src_v2_fname = field + "_src_v2_tmp.fits"
    src_t.write(src_v2_fname, overwrite=True)
else:
    #just use the SRC
    src_v2_fname=src
os.system('/soft/topcat/stilts tmatch2  join=all1 values1=ID values2=ID matcher=exact in1=%s in2=%s out=%s' % (mergeout,src_v2_fname,mergeout2))

t=Table.read(mergeout2)
finalname=mergeout2.replace('merged_src','final')
print 'Remove unnecessary columns'

cols=['RA_2','DEC_2','FLAG_OVERLAP_2','FLAG_CLEAN_2','id', 'ID_OPTICAL_2', 'ID_SPITZER_2','ID_2','EBV_2','FLAG_DEEP_2']
for fcol in t.colnames:
    if fcol.endswith("_fluxerr"):
        cols.append(fcol)
        cols.append(fcol[:-3])

for c in cols:
    if c in t.colnames:
        print c,
        sys.stdout.flush()
        t.remove_columns(c)

print
print 'Rename columns'
#t['Separation_1'].name='Separation'

for column in ['ID','RA','DEC','FLAG_OVERLAP','flag_clean', 'ID_OPTICAL', 'ID_SPITZER','EBV','FLAG_DEEP']:

    oldcol=column+'_1'
    if oldcol in t.colnames:
        t[oldcol].name=column
        print oldcol,'->',column,':',

# Read in XID+ catalogue and update the FIR data for sources that were re-run
xid_t = Table.read(xid)
in_fin = np.isin(t["Source_Name"], xid_t["Source_Name"])

# Hack for setting instrument name in flag_* to lower case - currently only needed for EN1
flags = [aa for aa in xid_t.colnames if aa.startswith("flag_")]
for flagc in flags:
    if "pacs" not in flagc.lower():
        xid_t[flagc].name = flagc.lower()
    else:
        print(flagc)
        flagc_split = flagc.split("_")
        flagc_new = flagc_split[0] + "_PACS_" + flagc_split[-1]
        xid_t[flagc].name = flagc_new

# Remove some duplicate columns
xid_cols = xid_t.colnames
spurious_cols = [aa for aa in xid_cols if aa.startswith("RA") or aa.startswith("Dec") or aa.startswith("Source_Name")]
for blah in spurious_cols:
    xid_cols.remove(blah)

# Update the XID+ columns into the final catalogue columns + any new columns
for col in xid_cols:
    if col in t.colnames:
        t[col][in_fin] = xid_t[col]
    else:
        # Make a new column
        t[col] = False
        t[col][in_fin] = xid_t[col]

print
print 'Remove whitespace padding:',
sys.stdout.flush()
stringcols=[n for (n,ty) in t.dtype.descr if 'S' in ty]
for c in stringcols:
    print c,
    sys.stdout.flush()
    t[c]=[s.rstrip() for s in t[c]]
print

print 'Updating the LR values and Position_from flag for subset that were later visually inspected'
# if field is not "bootes":
new_lr = Table.read(new_lr_path, character_as_bytes=False)
lr_vis = Table.read(lr_vis_path, format='ascii')

# Firstly, update the lr_th columns, setting a minimum value of 1.1*lr_th if needed?
raw_in_fin = np.isin(t["Source_Name"], new_lr["Source_Name"])
_, ind_t, ind_new_lr = np.intersect1d(np.array(t["Source_Name"]), np.array(new_lr["Source_Name"]), return_indices=True)

# Check that intersect1d is run correctly
assert np.sum(raw_in_fin) == len(ind_t), "Something has gone wrong in intersect1d!"

# Copy over the new LR values
t["lr_fin"][ind_t] = np.copy(new_lr["lr_fin"][ind_new_lr])

# Also set a minimum LR value? I don't think this is needed
# low_lr = t["lr_fin"] <= lr_th
# print("No. of low LR sources: {0}".np.sum(low_lr))
# t["lr_fin"][low_lr] = 1.1 * lr_th

# Update the Position_from for the subset that were visually inspected
vis_in_fin = np.isin(t["Source_Name"], lr_vis["Source_Name"])
t["Position_from"][vis_in_fin] = "Visual inspection"

print 'Remove -99s:',
dblcols=[n for (n,ty) in t.dtype.descr if ('f8' in ty or 'f4' in ty)]
for c in dblcols:
    print c,
    sys.stdout.flush()
    t[c]=np.where(t[c]==-99,np.nan,t[c])
print


print 'Sorting'
t.sort('RA')

print 'Finalize NoID values'
filter=(t['NoID']>0) & ~np.isnan(t['optRA'])
print 'Set',np.sum(filter),'sources with NoID set but some ID to NoID=0'
t['NoID']=np.where(filter,0,t['NoID'])
filter=~np.isnan(t['optRA']) & t['ID'].mask
print 'Set',np.sum(filter),'sources with optRA set but no ID to NoID=2'
t['NoID']=np.where(filter,2,t['NoID'])
t['optRA']=np.where(filter,np.nan,t['optRA'])
t['optDec']=np.where(filter,np.nan,t['optDec'])

print 'Writing to disk'
t.write(finalname,overwrite=True)
