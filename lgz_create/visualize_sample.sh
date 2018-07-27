#!/bin/bash

python $LGZPATH/utils/make_image_list.py $1
list=`echo $1 | sed -e 's/.fits/-list.txt/'`
echo $list
lines=`wc $list | awk '{print $1-1}'`
echo $lines
qsub -t 0-${lines} -v INFILE=$1,LGZPATH=$LGZPATH,IMAGEDIR=$IMAGEDIR $LGZPATH/lgz_create/visualize.qsub
