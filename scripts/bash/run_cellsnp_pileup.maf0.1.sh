#!/bin/bash

BAM=$1
BARCODES=$2
VCF=$3
OUTDIR=$4

singularity exec --bind /project/tishkofflab ~/tools/Demuxafy.sif cellsnp_pileup.py \
          -s $BAM \
          -b $BARCODES \
          -O $OUTDIR \
          -R $VCF \
          -p 20 \
          --minMAF 0.1 \
          --minCOUNT 20
