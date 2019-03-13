from astropy.io import fits
from astropy.wcs import WCS
import glob
import numpy as np

class ReadMaps(object):
    def __init__(self, gpath):
        self.files=glob.glob(gpath)
        self.hdulist=[]
        self.wcslist=[]
        self.shapelist=[]
        for f in self.files:
            hdu=fits.open(f)
            w=WCS(hdu[0].header)
            self.hdulist.append(hdu)
            self.wcslist.append(w)
            self.shapelist.append(hdu[0].data.shape)
        self.n=len(self.files)
        
    def find_pos(self,ra,dec,hdu=False):
        found=None
        dist=None
        fpos=None,None
        for i in range(self.n):
            w=self.wcslist[i]
            maxy,maxx=self.shapelist[i]
                                     
            centx=maxx/2
            centy=maxy/2
            x,y=w.wcs_world2pix(ra,dec,0)
            if x<0 or y<0 or x>=maxx or y>=maxy:
                continue
            cdist=np.sqrt((x-centx)**2.0+(y-centx)**2.0)
            if dist is None or cdist<dist:
                found=i
        if hdu and found is not None:
            return self.hdulist[found]
        else:
            return found
