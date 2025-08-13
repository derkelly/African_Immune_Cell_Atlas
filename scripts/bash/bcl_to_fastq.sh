#!/bin/bash

IN=$1
OUT=$2
SAMP=$3

module load bcl2fastq2/v2.20.0.422
~/bin/cellranger-7.1.0/cellranger mkfastq \
    --id=$OUT \
    --run=$IN \
    --csv=$SAMP
