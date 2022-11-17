from astropy.table import Table,vstack

for name in ['LGZ-multiple.fits','LGZ-comps.fits','LGZ-cat.fits']:

    t1=Table.read('Fall-rgzonly/'+name)
    t2=Table.read('lgz-0h/'+name)
    if name=='LGZ-cat.fits':
        t2['Hostbroken_prob'].name='Other_prob'
    t=vstack([t1,t2])
    t.write('Fall/'+name)
    
