#!/bin/bash
if [ $# -lt 3 ]
then
  echo "Usage: $0 bam pdf title"
  exit
fi

export t=`samtools view -H $1 | perl -e '$tot=0;while(<>){if (/@SQ/ && /LN:(\d+)/) {$tot+=$1}}print $tot;'`
bam2density --noheaders -p 0 --no_nh -s $1 | bamQC.R $2 $3 $t
