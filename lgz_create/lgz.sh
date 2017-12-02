#!/bin/bash

export PATH=/soft/Montage_v3.3/bin:$PATH
python $LGZPATH/make_overlays.py $INFILE $PBS_ARRAYID

