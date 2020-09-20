import os,sys,glob

from sample_rgz_output import run

root='/beegfs/lofar/mjh/rgz'
os.chdir(root)
g=glob.glob('P*')
for f in g:
    os.chdir(root+'/'+f)
    if not os.path.isfile('LGZ-cat-list.txt'):
        run('python $LGZPATH/utils/download_image_files_legacy.py LGZ-cat.fits')
    run('python ~/git/lotss-catalogue/lgz_create/make_overlays_legacy_inspect_rgzout.py LGZ-cat.fits 0 1000')
