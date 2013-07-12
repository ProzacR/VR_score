#!/bin/bash
#
# main fuction for training
# calls all scripts in proper order
#
# first need Final_table3.csv and structure files
# get it from create_set.sh and then prepare_set.sh
# bad part that it needs pach up with office (Kd values bit mess)
#
# VR


x=1
while [ -e result_${x}.csv ]; do
x=$((x + 1));
done
./train.sh > result_${x}.csv;
awk 'NR == 1' result_${x}.csv | sed -r 's/^.{4}//' > result_${x}c.csv;
awk '/^[0-9,.]+$/' result_${x}.csv >> result_${x}c.csv
R --no-save --args result_${x}c.csv < train.R
mv plot.png result_${x}.png
feh result_${x}.png &
