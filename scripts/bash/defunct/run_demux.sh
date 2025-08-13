#!/bin/bash

PILEUP=$1
VCF=$2
OUT=$3

conda activate muxlet

popscle demuxlet \
    --alpha 0 \
    --alpha 0.5 \
    --plp $PILEUP \
    --vcf <(zcat $VCF ) \
    --field GT \
    --out $OUT
