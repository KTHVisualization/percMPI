#!/bin/bash

dir=$1
awk '(NR==1)||(FNR>1)' $1*_timings_*.csv > $dir${dir%?}"_timings_all.txt" 

 
