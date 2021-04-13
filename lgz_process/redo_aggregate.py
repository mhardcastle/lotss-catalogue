import os,sys,glob

from sample_rgz_output import run

root='/beegfs/lofar/mjh/rgz'
os.chdir(root)
g=glob.glob('P*')
for f in g:
    os.chdir(root+'/'+f)
    ra,dec=f.replace('P','').split('+')
    ra=int(ra)
    dec=int(dec)
    field_l='13h' if (ra>110 and ra<280) else '0h'
    field_h='south' if (dec<33 or field_l=='0h') else 'north'
    run('python ~/git/lotss-catalogue/lgz_process/aggregate_lofgalzoo_dr2.py %s %s' % (field_l,field_h) )

    
