#!/bin/bash
set -e -o pipefail

#DEPENDENCIES samtools, freebayes

cd $1

mkdir -p freebayes_reference

echo "\
minimap2 -x map-hifi ../$2 filtered/all.fastq.gz -t 8 -a --secondary=no | samtools view -S -b -F 4 -F 0x800 > freebayes_reference/${1}.bam"
minimap2 -x map-hifi ../$2 filtered/all.fastq.gz -t 8 -a --secondary=no | samtools view -S -b -F 4 -F 0x800 > freebayes_reference/${1}.bam

echo "\
samtools sort freebayes_reference/${1}.bam -O BAM -o freebayes_reference/sorted_${1}.bam"
samtools sort freebayes_reference/${1}.bam -O BAM -o freebayes_reference/sorted_${1}.bam

echo "\
samtools index freebayes_reference/sorted_${1}.bam"
samtools index freebayes_reference/sorted_${1}.bam

echo "\
freebayes -f ../$2 freebayes_reference/sorted_${1}.bam -v freebayes_reference/$1.vcf"
freebayes -f ../$2 freebayes_reference/sorted_${1}.bam -v freebayes_reference/$1.vcf

