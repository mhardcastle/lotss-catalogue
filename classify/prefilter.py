import os
from astropy.table import Table
from subprocess import Popen
from time import sleep
import curses
import sys
import numpy as np

sys.stdout.write("\x1b]2;Classify\x07")
def displayfile(f):
#    return Popen('exec /usr/bin/display -resize 900x900 '+f+' >/dev/null 2>&1',shell=True)
    p=Popen('exec /usr/bin/eog '+f+' >/dev/null 2>&1',shell=True)
#    sleep(0.5)
#    cmd="wmctrl -r '"+os.path.basename(f)+"' -e 1,0,0,-1,-1"
#    print cmd
#    os.system(cmd)
    return p

def writeres():
    out=open(outfile,'w')
    for i,r in enumerate(results):
        out.write(r+','+t[i]['Source_Name']+'\n')
    out.close()


##### edit the following lines to choose the sample and possible options

dir='.'
sample='workflow'
options=('Send to LGZ', 'Accept ML match', 'No good match', 'Too zoomed in','Artefact','Uncatalogued host','Blend')

##### don't change anything else

wd=os.getcwd()
outfile=wd+'/'+sample+'.txt'
if os.path.isfile(outfile):
    results=[l.split(',')[0] for l in open(outfile).readlines()]
else:
    results=None
    
os.chdir(dir)
t=Table.read(sample+'.fits')
if results is None:
    results=['']*len(t)

assert(len(results)==len(t))

if len(sys.argv)>=3:
    start=int(sys.argv[1])
    end=int(sys.argv[2])
else:
    start=0
    end=len(t)
    
stdscr = curses.initscr()
curses.start_color()
curses.use_default_colors()
curses.cbreak()
for i in range(0, curses.COLORS):
    curses.init_pair(i, i, -1);

i=start
while results[i]!='':
    i+=1

while i<len(t) and i<end:
    r=t[i]
    stdscr.keypad(1)
    stdscr.clear()

    name=r['Source_Name']
    f=name+'_j.png'
    p=displayfile(f)
    sleep(0.8)
    os.system('wmctrl -a "Classify"')

    if results[i]=='':
        desc='Unclassified'
    else:
        desc=options[int(results[i])-1]
    instructions=[str(name),'LR '+str(r['lr_fin']),'(%i of %i)' % (i, len(t)),desc,'']
    instructions+=['(%i) %s' % (j+1,s) for j,s in enumerate(options)]
    instructions+=['','LEFT: back one','RIGHT: forward one', 'Q: quit','',"Your choice:"]
    attrs=[curses.color_pair(1)]*4+[None,]*(len(instructions)-4)
    if np.isnan(r['lr_fin']):
        attrs[6]=curses.A_REVERSE
    for j in range(len(instructions)):
        l=instructions[j]
        a=attrs[j]
        if a is not None:
            stdscr.addstr(4+j,2,l,a)
        else:
            stdscr.addstr(4+j,2,l)
    stdscr.refresh()

    valid=False
    quit=False
    while not(valid):
        key = stdscr.getch(4+j,2+len(instructions[-1])+1)

        valid=True
        if key==ord('q') or key==ord('Q'):
            quit=True
        elif key>ord('0') and key<=ord('9'):
            res=key-ord('0')
            results[i]=str(res)
            i+=1
        elif key==curses.KEY_LEFT:
            i-=1
        elif key==curses.KEY_RIGHT:
            i+=1
        else:
            valid=False

    writeres()
    p.terminate()
    if quit:
        break

curses.nocbreak()
stdscr.keypad(0)
curses.echo()

curses.endwin()
