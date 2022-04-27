from astropy.table import Table,vstack
import numpy as np

tm=Table.read('../mksp/LGZ_selection_all_radec.fits')
ts=Table.read('sources-v0.1.fits')
tc=Table.read('components-v0.1.fits')

ss=[]
origname=[]

for r in tm:
    compname=r['Source_Name'].rstrip()
    print compname
    filt=(tc['Component_Name']==compname)
    if np.sum(filt)==0:
        print compname,'not found!'
    else:
        cr=tc[filt][0]
        ps=cr['Parent_Source']
        pfilt=(ts['Source_Name']==ps)
        ss.append(ts[pfilt])
        origname.append(compname)

mkt=vstack(ss)
mkt['Original_Name']=origname

mkt.write('mksp-output.fits',overwrite=True)
