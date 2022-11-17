# bash initialization for catalogue code
# edit appropriately for your system

export LGZPATH=/home/mjh/git/lotss-catalogue
export IMAGEDIR=/data/lofar/DR2
export LOTSS_COMPONENT_CATALOGUE=/data/lofar/DR2/catalogues/LoTSS_DR2_v110_masked.srl.fits
#export LOTSS_COMPONENT_CATALOGUE=/data/lofar/DR2/catalogues/LoTSS_DR2_v100.srl.fits
export PATH=/soft/Montage_v3.3/bin:/soft/bin:$PATH

export PYTHONPATH=/soft/python/lib64/python2.7/site-packages:${LGZPATH}/utils:${LGZPATH}/dr2_catalogue:$PYTHONPATH
