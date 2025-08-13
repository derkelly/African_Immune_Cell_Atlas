#!/bin/bash

PLP=$1
VCF=$2
OUT=$3

singularity exec --bind /project/tishkofflab ~/tools/Demuxafy.sif popscle demuxlet \
    --plp $PLP \
    --vcf $VCF \
    --field GT \
    --out $OUT
