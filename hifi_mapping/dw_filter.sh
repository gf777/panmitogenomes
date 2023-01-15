#!/bin/bash
set -e -o pipefail

#DEPENDENCIES aws, pbindex, bam2fastq, minimap2, samtools, gzip, blast

PER=$3

cd $1

#VARIABLES
file=$(cat hifi_input.ls | sed -n ${SLURM_ARRAY_TASK_ID}p)
filename=$(basename -- "$file")

#DOWNLOAD
if [ ! -f raw/${filename} ]; then

    echo \
   "aws s3 cp $file raw"
    aws s3 cp $file raw
        
else

    echo "File raw/${filename} exists. Skipping."

fi

#CONVERT
echo "BAM>FASTQ for file raw/${filename}"

if [ ! -f fastq/${filename%.*}.fastq.gz ]; then
    
    pbindex raw/${filename}

    bam2fastq raw/${filename} -o fastq/${filename%.*}

else

    echo "File fastq/${filename%.*}.fastq.gz exists. Skipping..."

fi

#ALIGN
echo "Aligning fastq/${filename%.*}.fastq.gz..."

if [ ! -f alignments/${filename%.*}.bam ]; then

    echo "\
    minimap2 -x map-hifi ../$2 fastq/${filename%.*}.fastq.gz -t 8 -a --secondary=no | samtools view -S -b -F 4 -F 0x800 > alignments/${filename%.*}.bam"
    minimap2 -x map-hifi ../$2 fastq/${filename%.*}.fastq.gz -t 8 -a --secondary=no | samtools view -S -b -F 4 -F 0x800 > alignments/${filename%.*}.bam

else

    echo "File alignments/${filename%.*}.bam exists. Skipping..."

fi

#CONVERT
echo "BAM>FASTQ for file alignments/${filename%.*}.bam"

if [ ! -f alignments/${filename%.*}.fastq ]; then
    
    samtools fastq -F 0x900 -i alignments/${filename%.*}.bam > alignments/${filename%.*}.fastq
    
else

    echo "File alignments/${filename%.*}.fastq.gz exists. Skipping..."

fi

#FILTER (adapted from mitoVGP)
sed -n '1~4s/^@/>/p;2~4p' alignments/${filename%.*}.fastq > filtered/${filename%.*}.fasta

makeblastdb -in filtered/${filename%.*}.fasta -parse_seqids -dbtype nucl -out filtered/${filename%.*}.db

GSIZE=$(awk 'BEGIN {FS="\t"} $0 !~ ">" {sum+=length($0)} END {print sum}' ../${2%.*})

blastn -outfmt "6 sseqid slen qcovs" -query ../${2%.*} -db filtered/${filename%.*}.db | sort -k2 -nr | uniq | awk -v gsize="${GSIZE}" '{printf $0 "\t" gsize/$2*$3 "\n"}' | awk -v per="${PER}" '{if($4>per) printf $0 "\n"}' > filtered/filtered_${SLURM_ARRAY_TASK_ID}.out

printf "\n\nExtracting the following $(wc -l filtered/filtered_${SLURM_ARRAY_TASK_ID}.out | awk '{print $1}') reads by reference cover:\n\n"

sed -i 'Read ID\tRead length\tQry cover\tCorrected Qry cover' filtered/filtered_${SLURM_ARRAY_TASK_ID}.out

cat filtered/filtered_${SLURM_ARRAY_TASK_ID}.out | column -t

awk '{printf $1 "\n"}' filtered/filtered_${SLURM_ARRAY_TASK_ID}.out > filtered/filtered_${SLURM_ARRAY_TASK_ID}.ls

grep alignments/${filename%.*}.fastq -f filtered/filtered.ls -A3 --no-group-separator | gzip > filtered/${filename%.*}.fastq.gz

cat filtered/*fastq.gz > filtered/all.fastq.gz

#CLEANUP
gzip alignments/${filename%.*}.fastq
rm fastq/${filename%.*}.fastq.gz raw/*