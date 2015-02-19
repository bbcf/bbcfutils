#!/bin/sh

if [ $# -lt 4 ]
then
  echo "Usage: $0 [call|list] fasta [depth|bed] bam [samheader]"
  exit
fi

fasta=$2
depth=$3
bam=$4
if [ $# -gt 4 ]
then
    samhead=$5
fi
if [ $1 = "call" ]
then
    if [ -n "$samhead" ]
    then
        samtools reheader $samhead $bam | samtools mpileup -uDS -I -f $fasta - | bcftools view -vcg - | vcfutils.pl varFilter -D$depth
    else
        samtools mpileup -uDS -I -f $fasta $bam | bcftools view -vcg - | vcfutils.pl varFilter -D$depth
    fi
else
    if [ -n "$samhead" ]
    then
        samtools reheader $samhead $bam | samtools mpileup -uDS -I -f $fasta - | bcftools view -vcgl $depth - 
    else
        samtools mpileup -uDS -I -f $fasta $bam | bcftools view -vcgl $depth -
    fi
fi

