#!/bin/bash
set -e -o pipefail

#DEPENDENCIES samtools, freebayes

cd results/$1

mkdir freebayes

for bam in filtered/*
do

echo /"
samtools sort $bam -o filtered/sorted_$bam"
samtools sort $bam -o filtered/sorted_$bam

done

samtools merge filtered/sorted_* -o freebayes/merged_${1}.bam

samtools index merged_${1}.bam

freebayes -f ../$2 merged_${1}.bam -v freebayes/$1.vcf

#CLEANUP
rm filtered/sorted_*

