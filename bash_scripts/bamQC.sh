#!/bin/bash
if [ $# -lt 4 ]
then
  echo "Usage: $0 bam pdf genome_size title"
  exit
fi

bam2density --noheaders -p 0 --no_nh -s $1 | bamQC.R $2 $3 $4
