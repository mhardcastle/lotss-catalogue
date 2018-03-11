#!/bin/bash

export PATH=/soft/Montage_v3.3/bin:$PATH
export PYTHONPATH=$LGZPATH/utils:$PYTHONPATH
python $LGZPATH/lgz_create/visualize_catalogued.py $INFILE $PBS_ARRAYID
