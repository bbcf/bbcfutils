#!/bin/bash
bamname=$1
plotname=$2
bam2density -q 1 -p 0 --no_nh -s $bamname | cut -f 4 | sort -n | uniq -c | bamQC.R $plotname
