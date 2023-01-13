#!/bin/bash
set -e -o pipefail

#DEPENDENCIES samtools, freebayes

cd $1

mkdir -p freebayes

cat filtered/*fastq.gz > filtered/all.fastq.gz

echo "\
minimap2 -x map-hifi ../$2 filtered/all.fastq.gz -t 8 -a --secondary=no | samtools view -S -b -F 4 -F 0x800 > freebayes/${1}.bam"
minimap2 -x map-hifi ../$2 filtered/all.fastq.gz -t 8 -a --secondary=no | samtools view -S -b -F 4 -F 0x800 > freebayes/${1}.bam


echo "\
samtools sort freebayes/${1}.bam -O BAM -o filtered/sorted_${1}.bam"
samtools sort freebayes/${1}.bam -O BAM -o filtered/sorted_${1}.bam

echo "\
samtools merge filtered/sorted_${1}.bam -o freebayes/merged_${1}.bam"
samtools merge filtered/sorted_${1}.bam -o freebayes/merged_${1}.bam

echo "\
samtools index merged_${1}.bam"
samtools index merged_${1}.bam

echo "\
freebayes -f ../$2 merged_${1}.bam -v freebayes/$1.vcf"
freebayes -f ../$2 merged_${1}.bam -v freebayes/$1.vcf

