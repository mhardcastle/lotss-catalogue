from astropy.table import Table
from numpy import cos,pi
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import os
imsize=866
plot=False
import glob

# Modified to read the ellipse data from directories below the current one

#lines=[l.rstrip() for l in open('erik.txt').readlines()]

files=glob.glob('*/*ellipses.fits')
dirs=[]
sources=[]
for f in files:
    bits=f.split('/')
    dirs.append(bits[0])
    sources.append(bits[1].replace('-ellipses.fits',''))

files=[]
for sourcename,directory in zip(sources,dirs):
    print 'Doing %s in directory %s' % (sourcename, directory)
    os.chdir('/data/lofar/DR2/RGZ/'+directory)
    manifest=open(sourcename+'-manifest.txt').readlines()[0]
    bits=manifest.rstrip().split(',')
    ra=float(bits[5])
    dec=float(bits[6])
    size=float(bits[7])
    scale=imsize/size # pixels per arcsec
    png=bits[1]
    if plot:
        img=plt.imread(png)
        plt.imshow(img)
        ax=plt.gca()
    # now do ellipses
    with open(sourcename+'-ellipses.txt','w') as outfile:
        t=Table.read(sourcename+'-ellipses.fits')
        for r in t:
            x=imsize/2+scale*3600*(ra-r['RA'])*cos(dec*pi/180)
            y=imsize/2+scale*3600*(dec-r['DEC'])
            if x<0 or y<0 or y>imsize or x>imsize:
                continue
            outfile.write('%f %f %f %f %f\n' % (x,y,r['Maj']*scale,r['Min']*scale,360-r['PA']))
            if plot:
                plt.scatter(x,y,marker='o',color='green',s=100)
                ellipse=Ellipse((x,y),height=2*r['Maj']*scale,width=2*r['Min']*scale,angle=360-r['PA'])
                ax.add_artist(ellipse)
                ellipse.set_facecolor('green')
                ellipse.set_alpha(0.5)
    if plot: plt.show()
    files=files+[directory+'/'+f for f in [bits[1],bits[2],bits[3],sourcename+'-ellipses.txt']]
        
print files
os.chdir('/data/lofar/DR2/RGZ/')
os.system('tar cvf erik.tar '+' '.join(files))
