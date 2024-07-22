#!/bin/sh
i=$1
j=${i%.*}
CAT contigs -c $i -o $j'.cat' --force \
-d /public3/home/sc52870/soft/CAT_prepare_20210107/2021-01-07_CAT_database/ \
-t /public3/home/sc52870/soft/CAT_prepare_20210107/2021-01-07_taxonomy/ \
-n 64

CAT add_names -i $j'.cat.contig2classification.txt' \
-o $j'.contig2classification.addname.txt' \
-t /public3/home/sc52870/soft/CAT_prepare_20210107/2021-01-07_taxonomy/

rm $j'.cat.'*
