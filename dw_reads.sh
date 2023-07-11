#!/bin/bash

while read -r sample others; do

mkdir -p results/${sample}

sbatch -pvgl -c8 --output=dw_reads.%A_%a.out --error=dw_reads.%A_%a.out aws s3 --no-sign-request cp \
    --recursive \
    s3://human-pangenomics/submissions/1E2DD570-3B26-418B-B50F-5417F64C5679--HIFI_DEEPCONSENSUS/${sample}/raw_data/PacBio_HiFi/DeepConsensus/ \
    results/${sample}

done<$1
