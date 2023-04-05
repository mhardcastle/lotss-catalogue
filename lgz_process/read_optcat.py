# test code

from __future__ import print_function
from astropy.table import Table
from multiprocessing import Pool
from tqdm import tqdm
import glob

g=glob.glob('optical_hpix/*.fits')

hp=[int(f.replace('optical_hpix/','').replace('.fits','')) for f in g]

def load(h):
    return Table.read('optical_hpix/%i.fits' % h)

p=Pool(16)
od={}
for i,result in enumerate(tqdm(p.imap(load,hp),total=len(hp))):
    od[hp[i]]=result
