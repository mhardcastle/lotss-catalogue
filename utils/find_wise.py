from __future__ import print_function
import os
import glob
from astropy.io import fits
from astropy.wcs import WCS

class WISE(object):
    def __init__(self):

        imagedir=os.environ['IMAGEDIR']

        g=glob.glob(imagedir+'/downloads/*w1*fits')
        self.wd={}
        self.wcs={}
        for f in g:
            name=f.split('/')[-1]
            header=fits.open(f)[0].header
            self.wd[name]=header

        for k in self.wd:
            self.wcs[k]=WCS(self.wd[k])

    def find_pos(self,ra,dec):
        found=[]
        for f in self.wd:
            x,y=self.wcs[f].wcs_world2pix(ra,dec,0)
            if x>=0 and y>=0 and x<=self.wd[f]['NAXIS1'] and y<=self.wd[f]['NAXIS2']:
                found.append((f,x,y))
        if len(found)==0:
            return None
        elif len(found)==1:
            return found[0][0]
        else:
            # don't know if overlapping images exist but do the right thing
            dists=[]
            maxdist=None
            selected=None
            for f,x,y in found:
                dist=(x-self.wd[f]['NAXIS1']/2)**2 + (y-self.wd[f]['NAXIS2'])**2
                if maxdist is None or dist<maxdist:
                    maxdist=dist
                    selected=f
            return selected

if __name__=='__main__':
    
    w=WISE()
    print(w.find_pos(3.84002098832,26.215683769))
    
