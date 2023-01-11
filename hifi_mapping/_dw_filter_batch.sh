#!/bin/bash
set -e -o pipefail

# $1 = metadata, $2 = reference, $3 = filter threshold (% breadth of coverage)

#WD
mkdir -p results
ln -f $2 results
cd results

#INDEX REFERENCE
reference=$(basename -- "$2")
if [ ! -f $reference.mmi ]; then

    echo "Indexing $reference"

    minimap2 -d $reference.mmi $reference

fi

#SBATCH SAMPLE JOBS
{
    read
    while read p; do

        IFS=',' read -r -a LINE <<< "$p"      
        ID=${LINE[0]}
        echo "processing sample $ID"

        mkdir -p $ID/log $ID/raw $ID/fastq $ID/alignments $ID/filtered

        #COLLECT READ PATHS
        IFS=';' read -r -a BAM <<< "${LINE[5]}"

        rm -f $ID/hifi_input.ls

        for i in ${BAM[@]}; do
            echo $i
            echo $i >> $ID/hifi_input.ls
        done

        FILE_NUM=$(wc -l < $ID/hifi_input.ls)

        echo "\
        sbatch -pvgl --array=1-${FILE_NUM} --output=$ID/log/dw.%A_%a.out --error=$ID/log/dw.%A_%a.out --cpus-per-task=8 ../dw_filter.sh $ID $reference.mmi $3"
        sbatch -pvgl --array=1-${FILE_NUM} --output=$ID/log/dw.%A_%a.out --error=$ID/log/dw.%A_%a.out --cpus-per-task=8 ../dw_filter.sh $ID $reference.mmi $3 | awk '{print $4}' >> dw.jid

    done
}<../$1