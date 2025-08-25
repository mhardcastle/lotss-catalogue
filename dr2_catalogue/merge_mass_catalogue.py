# Merge Ken's galaxy mass catalogue

from __future__ import print_function
import glob
from astropy.table import Table,hstack,MaskedColumn,join
import numpy as np
from astropy.io import fits

def add_descs(tj,outname):
    
    print('Add descriptions')
    lines=open('README-mass.md').readlines()
    for l in lines:
        l=l.replace('\\','')
        if '|' not in l:
            continue
        if 'Column' in l:
            continue
        if 'Flag value' in l:
            break
        bits=l.split('|')
        if len(bits)>1:
            c=bits[1].lstrip().rstrip()
            u=bits[2].lstrip().rstrip()
            d=bits[3].lstrip().rstrip()
            if c[0]==':': continue

            if c in tj.colnames:
                tj[c].description=d
                tj[c].units=u

    print('Writing')
    tj.write(outname,overwrite=True)

    print('Fix the descriptions again')
    hdu=fits.open(outname)

    nf=hdu[1].header['TFIELDS']
    for i in range(1,nf+1):
        c=hdu[1].header['TTYPE%i' %i]
        try:
            hdu[1].header['TCOMM%i' %i]=tj[c].description
            hdu[1].header['TUNIT%i' %i]=tj[c].units
        except AttributeError:
            print(c,'has no units or description')

    hdu.writeto(outname,overwrite=True)

if __name__=='__main__':

    release=sorted(glob.glob('/beegfs/lofar/mjh/rgz/combined-release-v*-LM.fits'))[-1]

    kdcat='/beegfs/lofar/duncan/combined_release_v1.2_opt_mass_cols.fits'
    #kdcat='/beegfs/lofar/mjh/rgz/fixed_mass_table.fits'

    print('Loading...')
    print('...',release)
    t=Table.read(release)
    print('...',kdcat)
    tkd=Table.read(kdcat)


    # Remove duplicate columns in original table
    print('Deduping')

    for k in tkd.colnames[1:]: # not Source_Name
        if k in t.colnames:
            del(t[k])
        if k.upper() in t.colnames:
            del(t[k.upper()])

    # Remove luptitudes

    print('Deluptituding')
    for k in t.colnames:
        if 'lupt' in k:
            del(t[k])

    print('Joining')
    tj=Table(join(t,tkd,join_type='left'),masked=True)

    tj['flag_mass_bool']=np.where((tj['flag_mass']==1) & ~tj['flag_mass'].mask,True,False)
    #tj['flag_mass_bool'].mask=tj['flag_mass'].mask

    del(tj['flag_mass'])
    tj['flag_mass_bool'].name='flag_mass'

    outname=release.replace('.fits','_opt_mass.fits')
    add_descs(tj,outname)
    
