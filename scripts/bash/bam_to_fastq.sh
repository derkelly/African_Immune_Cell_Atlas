#!/bin/bash

BAM=$1
OUTDIR=$2

~/bin/cellranger-7.1.0/lib/bin/bamtofastq \
    --nthreads=8 \
    $BAM \
    $OUTDIR
