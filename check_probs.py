from astropy.table import Table
import os.path

t=Table.read('drop_noprobs.fits')
compt=Table.read('HETDEX-LGZ-comps-v0.3.fits')

for r in t:
    name=r['Source_Name']
    print name
    regionfile=name+'.png'
    assert(os.path.isfile(regionfile))
    old_com=compt['Source_Name']==r['Source_Name']
    comps=compt[old_com]
    for r2 in comps:
        compfile='/data/lofar/mjh/hetdex_v3/lgz/subset_with_id_dir/'+r2['Comp_Name']+'_PS.png'
        print compfile
        if os.path.isfile(compfile):
            break
    else:
        print 'Argh argh argh'
        continue
    outfile=name+'_join.png'
    print regionfile,compfile
    os.system('convert +append '+regionfile+' '+compfile+' '+outfile)
    
