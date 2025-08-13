#!/bin/bash

PLP=$1
NSAMP=$2
OUT=$3

singularity exec --bind /project/tishkofflab ~/tools/Demuxafy.sif popscle freemuxlet \
    --plp $PLP \
    --out $OUT \
    --nsample $NSAMP
