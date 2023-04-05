from __future__ import print_function
import os
from astropy.table import Table
from subprocess import Popen
from time import sleep
import curses
import sys
import numpy as np
import MySQLdb as mdb
import MySQLdb.cursors as mdbcursors
import matplotlib.pyplot as plt 
import matplotlib.image as img 
import glob

sys.stdout.write("\x1b]2;Classify\x07")

def get_next():
    cur.execute('lock table %s write' % table)
    cur.execute('select * from %s where user is NULL order by object limit 1' % table )
    res=cur.fetchall()
    if len(res)==0:
        object=None
    else:
        object=res[0]['object']
        cur.execute('update %s set user="%s" where object="%s"' % (table,user,object))
    cur.execute('unlock tables')
    return res[0]
                    
    
##### edit the following lines to choose the sample and possible options

table='blend_filter_Spring'
dir='/beegfs/lofar/mjh/rgz/Spring/ngp_prefilter/'  #%s-prefilter' % table
os.chdir(dir)

g=glob.glob('*.fits')
assert(len(g)==1)
sample=g[0]

options=('Drop all IDs', 'Accept main ID and discard others', 'Accept main ID but split off Gaussians with plausible ID', 'Split up into all Gaussians with IDs','Send to blend','Artefact')
user=os.getenv('USER')

##### don't change anything else

con=mdb.connect('127.0.0.1', 'prefilter_user', 'WQ98xePI', 'prefilter', cursorclass=mdbcursors.DictCursor, autocommit=True)
cur = con.cursor()

plt.figure(figsize=(10,5),facecolor='black')

# get my old results if any

cur.execute('select * from %s where user="%s"' % (table, user))

results=list(cur.fetchall())

#print results
i=len(results)
    
t=Table.read(sample)
if 'ra' in t.colnames:
    t['ra'].name='optRA'
    t['dec'].name='optDec'

stdscr = curses.initscr()
curses.start_color()
curses.use_default_colors()
curses.cbreak()
for i in range(0, curses.COLORS):
    curses.init_pair(i, i, -1);

i=len(results)
while True:
    if i<len(results):
        # we are looking back at an old one
        name=results[i]['object']
        classification=results[i]['classification']
        gaussian_ids=results[i]['gaussian_ids']
    else:
        # need a new target
        result=get_next()
        name=result['object']
        if name is None:
            break
        gaussian_ids=result['gaussian_ids']
        classification=None
    r=t[t['Source_Name']==name][0]
    stdscr.keypad(1)
    stdscr.timeout(10)
    stdscr.clear()

    f=name+'_j.png'
    im=img.imread(f)
    plt.clf()
    plt.imshow(im)
    plt.axis('off')
    plt.tight_layout()
    plt.ion()
    plt.show()
    
    cur.execute('select count(object) from %s where classification is NULL' % table)
    res=cur.fetchall()
    remaining=res[0]['count(object)']
    plt.pause(0.1)

    #sleep(1.5)
    #os.system('wmctrl -a "Classify"')

    if classification is None:
        desc='Unclassified'
    else:
        desc=options[classification-1]
    instructions=[str(name),'(%i done, %i to do)' % (i, remaining),desc,'']
    instructions+=['(%i) %s' % (j+1,s) for j,s in enumerate(options)]
    instructions+=['','LEFT: back one','RIGHT: forward one', 'Q: quit','',"Your choice:"]
    attrs=[curses.color_pair(1)]*4+[curses.A_BOLD]*len(options)+[None,]*6
    lrvalid=True
    if np.isnan(r['optRA']):
        attrs[5]=None
        attrs[6]=None
        lrvalid=False
    if gaussian_ids==0:
        attrs[7]=None
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
            if (res==2 or res==3) and lrvalid==False:
                valid=False
                continue
            if (res==4 and gaussian_ids==0):
                valid=False
                continue
            if res>len(options):
                valid=False
                continue
            if i==len(results):
                results.append({'object':name,'classification':res,'gaussian_ids':gaussian_ids})
            else:
                results[i]['classification']=res
            cur.execute('update %s set classification=%i where object="%s"' % (table,res,name))
            i+=1
        elif key==curses.KEY_LEFT and i>0:
            if i==len(results):
                results.append({'object':name,'classification':None,'gaussian_ids':gaussian_ids})
            i-=1
        elif key==curses.KEY_RIGHT:
            if i==len(results):
                results.append({'object':name,'classification':None,'gaussian_ids':gaussian_ids})
            i+=1
        else:
            valid=False
            plt.pause(0.001)

    #p.terminate()
    if quit:
        break

curses.nocbreak()
stdscr.keypad(0)
curses.echo()

curses.endwin()

# hand back any objects we didn't classify
cur.execute('update %s set user=NULL where user="%s" and classification is NULL' % (table,user))

con.close()
