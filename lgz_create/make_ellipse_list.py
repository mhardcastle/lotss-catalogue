from astropy.table import Table
from numpy import cos,pi
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import os
imsize=866
plot=False

lines=[l.rstrip() for l in open('erik.txt').readlines()]

files=[]
for l in lines:
    if l:
        bits=l.split()
        sourcename=bits[1].replace('_PS','')
        print 'Doing %s (%s)' % (sourcename,bits[0][:-1])
        manifest=open(sourcename+'-manifest.txt').readlines()[0]
        bits=manifest.rstrip().split(',')
        ra=float(bits[4])
        dec=float(bits[5])
        size=float(bits[6])
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
        files=files+[bits[1],bits[2],sourcename+'-ellipses.txt']
        
print files
os.system('tar cvf hugh.tar '+' '.join(files))
