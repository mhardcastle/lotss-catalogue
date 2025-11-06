import glob
import sys
import os
from subprocess import check_output

queued=[]
queue=check_output('qstat -a',shell=True,universal_newlines=True).split('\n')
for l in queue:
    if 'hpmap-' in l:
        bits=l.split()
        field=bits[3][6:]
        print('Found',field,'already in queue')
        queued.append(field)

os.chdir('/beegfs/lofar/DR3/healpix_mosaics')
g=glob.glob('*')

for f in g:
    if os.path.isdir(f) and not os.path.isfile(f+'/hptable.fits') and os.path.isfile(f+'/mosaic-blanked.fits') and f not in queued:
        command=f'qsub -q core32 -N hpmap-{f} -v MOSAIC={f} /home/mjh/git/lotss-catalogue/dr3_tests/mos2hpmap_hp.qsub'
        print(command)
