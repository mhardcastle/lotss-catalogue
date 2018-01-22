import os
import numpy as np
import glob

def find_bbox(t):
    # given a table t find the bounding box of the ellipses for the regions

    boxes=[]
    for r in t:
        a=r['Maj']/3600.0
        b=r['Min']/3600.0
        th=(r['PA']+90)*np.pi/180.0
        dx=np.sqrt((a*np.cos(th))**2.0+(b*np.sin(th))**2.0)
        dy=np.sqrt((a*np.sin(th))**2.0+(b*np.cos(th))**2.0)
        boxes.append([r['RA']-dx/np.cos(r['DEC']*np.pi/180.0),
                      r['RA']+dx/np.cos(r['DEC']*np.pi/180.0),
                      r['DEC']-dy, r['DEC']+dy])

    boxes=np.array(boxes)

    minra=np.nanmin(boxes[:,0])
    maxra=np.nanmax(boxes[:,1])
    mindec=np.nanmin(boxes[:,2])
    maxdec=np.nanmax(boxes[:,3])
    
    ra=np.mean((minra,maxra))
    dec=np.mean((mindec,maxdec))
    size=1.2*3600.0*np.max((maxdec-mindec,(maxra-minra)*np.cos(dec*np.pi/180.0)))
    return ra,dec,size

def get_mosaic_name(name):
    globst=os.environ['IMAGEDIR']+'/mosaics/'+name.rstrip()+'*'
    g=glob.glob(globst)
    if len(g)==1:
        return g[0]
    elif len(g)==0:
        raise RuntimeError('No mosaic called '+name)
    else:
        raise RuntimeError('Mosaic name ambiguous')

