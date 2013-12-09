#!/bin/sh

if [ $# -lt 3 ]
then
  echo "Usage: `` fasta depth bam [samheader]"
  exit
fi

fasta=$1
depth=$2
bam=$3
if [ $# -gt 3 ]
then
  samtools reheader $4 $bam | samtools mpileup -uDS -I -f $fasta - | bcftools view -vcg - | vcfutils.pl varFilter -D$depth
else
  samtools mpileup -uDS -I -f $fasta $bam | bcftools view -vcg - | vcfutils.pl varFilter -D$depth
fi

