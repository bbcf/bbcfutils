#!/bin/sh

BRICKS=$1
frags=$2".tmp"
resfile=$3

sed 1d $2 > $frags

## return only frags with BRICKS (those correspond to frags located on the other chromosomes)
## (for now the found BRICKS are in the viewpoint's chromosome only)
sed 1d $BRICKS | awk '{print $1"\t"$2"\t"$3"\t"(-log($5)/log(10))}'|intersectBed -a $frags -b stdin -wao | awk '{if($5~/chr/){print $1"\t"$2"\t"$3"\t"$8}}'> $resfile
