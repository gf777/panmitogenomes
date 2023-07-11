#!/bin/bash
set -e -o pipefail

#DEPENDENCIES gfastats, samtools, freebayes

cd $1

mkdir -p freebayes_self

echo "\
gfastats ../../${2} $(grep $1 ../../${2} | tr -d ">") -o $1.fasta"
gfastats ../../${2} $(grep $1 ../../${2} | tr -d ">") -o $1.fasta

echo "\
minimap2 -d $1.fasta.mmi $1.fasta"
minimap2 -d $1.fasta.mmi $1.fasta

echo "\
minimap2 -x map-hifi $1.fasta.mmi filtered/*.fastq.gz -t 8 -a --secondary=no | samtools view -S -b -F 4 -F 0x800 > freebayes_self/${1}.bam"
minimap2 -x map-hifi $1.fasta.mmi filtered/*.fastq.gz -t 8 -a --secondary=no | samtools view -S -b -F 4 -F 0x800 > freebayes_self/${1}.bam

echo "\
samtools sort freebayes_self/${1}.bam -O BAM -o freebayes_self/sorted_${1}.bam"
samtools sort freebayes_self/${1}.bam -O BAM -o freebayes_self/sorted_${1}.bam

echo "\
samtools index freebayes_self/sorted_${1}.bam"
samtools index freebayes_self/sorted_${1}.bam

echo "\
freebayes -f $1.fasta freebayes_self/sorted_${1}.bam -v freebayes_self/$1.vcf --min-alternate-fraction 0.3"
freebayes -f $1.fasta freebayes_self/sorted_${1}.bam -v freebayes_self/$1.vcf --min-alternate-fraction 0.3

