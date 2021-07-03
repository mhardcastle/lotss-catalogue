import matplotlib.pyplot as plt 
import matplotlib.image as img 
import glob
from time import sleep

g=glob.glob('*_j.png')

plt.figure(figsize=(10,5),facecolor='black')
for f in g:
    im=img.imread(f)
    plt.imshow(im)
    plt.axis('off')
    plt.tight_layout()
    plt.show(block=False)
    while True:
        sleep(0.1)
        plt.pause(0.001)
    plt.clf()
    
