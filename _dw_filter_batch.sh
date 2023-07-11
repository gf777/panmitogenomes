#!/bin/bash
set -e -o pipefail

# $1 = metadata, $2 = reference, $3 = filter threshold (% breadth of coverage), $4 = individual assembly files

#WD
mkdir -p results
ln -f $2 results
ln -f $4 results
cd results

#INDEX REFERENCE
reference=$(basename -- "$2")
if [ ! -f $reference.mmi ]; then

    echo "Indexing $reference"

    minimap2 -d $reference.mmi $reference

fi

#SBATCH SAMPLE JOBS
{
#    read
    while read p; do

        IFS=',' read -r -a LINE <<< "$p"      
        ID=${LINE[0]}
        echo "processing sample $ID"

        mkdir -p $ID/log $ID/raw $ID/fastq $ID/alignments $ID/filtered

        #COLLECT READ PATHS
#         IFS=';' read -r -a BAM <<< "${LINE[1]}"
# 
#         rm -f $ID/hifi_input.ls
# 
#         for i in ${BAM[@]}; do
#             echo $i
#             echo $i >> $ID/hifi_input.ls
#         done

        ls ${LINE[1]} > hifi_input.ls

        FILE_NUM=$(wc -l < $ID/hifi_input.ls)

        echo "\
        sbatch -pvgl --array=1-${FILE_NUM} --output=$ID/log/dw.%A_%a.out --error=$ID/log/dw.%A_%a.out --cpus-per-task=32 ../dw_filter.sh $ID $reference.mmi $3"
        sbatch -pvgl --array=1-${FILE_NUM} --output=$ID/log/dw.%A_%a.out --error=$ID/log/dw.%A_%a.out --cpus-per-task=32 ../dw_filter.sh $ID $reference.mmi $3 | awk '{print $4}' > ${ID}/dw_filter.jid
        
        echo "\
        sbatch -pvgl --dependency=afterok:$(cat ${ID}/dw_filter.jid) --output=$ID/log/freebayes_reference.%A_%a.out --error=$ID/log/freebayes_reference.%A_%a.out --cpus-per-task=8 ../freebayes_reference.sh $ID $reference"
        sbatch -pvgl --dependency=afterok:$(cat ${ID}/dw_filter.jid) --output=$ID/log/freebayes_reference.%A_%a.out --error=$ID/log/freebayes_reference.%A_%a.out --cpus-per-task=8 ../freebayes_reference.sh $ID $reference | awk '{print $4}' > ${ID}/freebayes_reference.jid  
        
        echo "\
        sbatch -pvgl --dependency=afterok:$(cat ${ID}/dw_filter.jid) --output=$ID/log/freebayes_self.%A_%a.out --error=$ID/log/freebayes_self.%A_%a.out --cpus-per-task=8 ../freebayes_self.sh $ID $4"
        sbatch -pvgl --dependency=afterok:$(cat ${ID}/dw_filter.jid) --output=$ID/log/freebayes_self.%A_%a.out --error=$ID/log/freebayes_self.%A_%a.out --cpus-per-task=8 ../freebayes_self.sh $ID $4     

    done
}<../$1