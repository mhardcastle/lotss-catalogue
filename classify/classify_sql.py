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
    cur.execute('select object from %s where user is NULL order by object limit 1' % table )
    res=cur.fetchall()
    if len(res)==0:
        object=None
    else:
        object=res[0]['object']
        cur.execute('update %s set user="%s" where object="%s"' % (table,user,object))
    cur.execute('unlock tables')
    return object
                    
    
##### edit the following lines to choose the sample and possible options

dir='/beegfs/lofar/mjh/classify-2arcmin'
bits=dir.split('/')
table=bits[-1].replace('-','_')
os.chdir(dir)

g=glob.glob('sources*.fits')
assert(len(g)==1)
sample=g[0]
user=os.getenv('USER')

classifications=(('Basic',
                  (('fri','FRI'),
                   ('frii','FRII'),
                   ('headtail','Head-tail'),
                   ('nat','Narrow-angle tail (NAT)'),
                   ('wat','Wide-angle tail (WAT)'),
                   ('relaxed','Relaxed double'))),
                  ('Features',
                   (('bent','Bent'),
                    ('wings','Wings'),
                    ('xshaped','X-shaped'),
                    ('onesided','One-sided'),
                    ('multiple','Multiple hotspots'),
                    ('continuous','Continuous jet'),
                    ('restart','Restart'))),
                   ('Special',
                    (('cluster','Cluster'),
                     ('fluff','Fluff'),
                     ('revisit','Revisit'))))
                    
##### don't change anything else

keys='0123456789abcdeghijklmnoprsuvwxyzABCDEGHIJKLMNOPRSUVWXYZ'
labels=[]
names=[]
letters={}
classes=[]
maxlen=0
i=0
for clss,values in classifications:
    classes.append(clss)
    for label,name in values:
        labels.append(label)
        names.append(name)
        if len(name)>maxlen:
            maxlen=len(name)
        letters[label]=keys[i]
        i+=1

if len(sys.argv)>1:
    scommand=sys.argv[1]
else:
    scommand=None

t=Table.read(sample)

con=mdb.connect('127.0.0.1', 'classify_user', 'UI98rbYH', 'classify', cursorclass=mdbcursors.DictCursor, autocommit=True)
cur = con.cursor()

if scommand=='create':
    print('Creating table!')
    command='DROP TABLE '+table
    try:
        cur.execute(command)
    except mdb.Error as e:
        print(e)
        pass
    command='CREATE TABLE '+table+' (id int not null auto_increment primary key, user varchar(16), source varchar(32), done int, '
    for l in labels:
        command+=l+' int, '
    command=command[:-2]+')'
    print(command)
    cur.execute(command)

command='select * from %s where user="%s"' % (table,user)
cur.execute(command)
res=cur.fetchall()
if len(res)==0:
    scommand='populate' # force populate if table exists but has no entries for user'
    
if scommand=='create' or scommand=='populate':
    print('Populating table!')
    command='DELETE FROM '+table+' WHERE user="'+user+'"'
    cur.execute(command)
    for r in t:
        command='INSERT INTO '+table+' VALUES (NULL, "%s", "%s", 0, ' % (user,r['Source_Name'])
        for l in labels:
            command+='0, '
        command=command[:-2]+')'
        cur.execute(command)

plt.figure(figsize=(8,8),facecolor='black')

command='select source from %s where user="%s" and done=0' % (table,user)
cur.execute(command)
notdone=[r['source'] for r in cur.fetchall()]

for i in range(len(t)):
    name=t[i]['Source_Name']
    if name in notdone:
        break
else:
    raise RuntimeError('No sources need doing')

    
stdscr = curses.initscr()
curses.start_color()
curses.use_default_colors()
curses.cbreak()
for j in range(0, curses.COLORS):
    curses.init_pair(j, j, -1);

show_ellipses=True
while True:
    start=2
    cstart=2
    r=t[i]
    name=r['Source_Name']
    cur.execute('select * from %s where user="%s" and source="%s"' % (table,user,name))
    res=cur.fetchall()
    db=res[0] # a dictionary of key/value pairs for this source
    stdscr.keypad(1)
    stdscr.timeout(10)
    stdscr.clear()

    f=name+('_C.png' if show_ellipses else '_Cp.png')
    im=img.imread(f)
    plt.clf()
    plt.imshow(im)
    plt.axis('off')
    plt.tight_layout()
    plt.ion()
    plt.show()
    
    plt.pause(0.1)

    instructions=[str(name)+' (%i/%i)' % (i,len(t)),'Origin '+r['Created']]
    for j in range(len(instructions)):
        stdscr.addstr(start+j,cstart,instructions[j],curses.color_pair(1))

    start+=j+2 # starting row for classifications
    colwidth=maxlen + 5
    colsavail=curses.COLS//colwidth
    rowsavail=curses.LINES-start-7
    rowsneeded=1+((len(names)+len(classes))//colsavail)
    if rowsavail < rowsneeded:
        curses.nocbreak()
        stdscr.keypad(0)
        curses.echo()
        curses.endwin()
        raise RuntimeError('Not enough terminal real estate')
    row=start
    col=0
    for clss,values in classifications:
        stdscr.addstr(row,cstart+col*colwidth,'=== '+clss+' ===')
        row+=1
        if (row-start)>rowsneeded:
            row=start
            col+=1
        for label,vname in values:
            key=letters[label]
            if db[label]==1:
                stdscr.addstr(row,cstart+col*colwidth,'(%s) %s' % (key,vname),curses.A_REVERSE)
            else:
                stdscr.addstr(row,cstart+col*colwidth,'(%s) %s' % (key,vname))
            row+=1
            if (row-start)>rowsneeded:
                row=start
                col+=1
    
    start+=rowsneeded+1
                
    instructions=['','LEFT: back one','RIGHT or CR: forward one', 'F: see FITS','T: toggle labels','Q: quit','',"Your choice:"]
    for j in range(len(instructions)):
        l=instructions[j]
        stdscr.addstr(start+j,2,l)
    stdscr.refresh()
    valid=False
    quit=False
    newi=i
    while not(valid):
        key = stdscr.getch(start+j,2+len(instructions[-1])+1)
        valid=True
        if key==ord('q') or key==ord('Q'):
            quit=True
        elif key==ord('T') or key==ord('t'):
            show_ellipses=not(show_ellipses)
        elif key==ord('f') or key==ord('F'):
            os.system('/soft/bin/ds9 '+name+'.fits &')
        elif key==curses.KEY_LEFT and i>0:
            newi=i-1
        elif key==curses.KEY_RIGHT or key==curses.KEY_ENTER or key==13 or key==10:
            newi+=1
        elif key>=ord('0') and key<=ord('z'):
            # check if this is one of the keys in the section
            for k in letters:
                dbkey=letters[k]
                if key==ord(dbkey):
                    break
            else:
                valid=False
                continue
            newvalue=(0 if db[k]==1 else 1)
            command='UPDATE %s set %s=%i where user="%s" and source="%s"' % (table,k,newvalue,user,name)
            cur.execute(command)
        else:
            valid=False

        plt.pause(0.001)
        if newi!=i or quit:
        # source is done
            command='UPDATE %s set done=1 WHERE user="%s" and source="%s"' % (table,user,name)
            cur.execute(command)
            i=newi

    if quit:
        break

curses.nocbreak()
stdscr.keypad(0)
curses.echo()

curses.endwin()

con.close()
