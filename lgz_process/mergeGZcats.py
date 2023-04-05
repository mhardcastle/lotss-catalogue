from astropy.table import Table,vstack

# changed 29/11/22 to pick up weighted catalogue

weighted=True

for name in ['LGZ-multiple.fits','LGZ-comps.fits','LGZ-cat.fits']:

    if weighted:
        t1=Table.read('Fall-rgzonly/weighted-'+name)
    else:
        t1=Table.read('Fall-rgzonly/'+name)
    t2=Table.read('lgz-0h/'+name)
    if name=='LGZ-cat.fits':
        t2['Hostbroken_prob'].name='Other_prob'
    t=vstack([t1,t2])
    t.write('Fall/'+name,overwrite=True)
    
