#!bin/bash

ID=$1
FASTQP=$2
SAMP=$3

~/bin/cellranger-7.1.0/cellranger count \
    --id=$ID \
    --transcriptome=/project/tishkofflab/data_public/reference/refdata-gex-GRCh38-2020-A \
    --fastqs=$FASTQP \
    --sample=$SAMP
