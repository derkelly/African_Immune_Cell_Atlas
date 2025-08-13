#!/bin/bash

IDS=$1
GENOS=$2
OUT=$3

bcftools view \
    --samples-file $IDS \
    -R /project/tishkofflab/projects/SingleCellRNA/sample_info/refdata-gex-GRCh38-2020-A.genes.merge.bed \
    -Ou \
    $GENOS | \
    bcftools +fill-tags -Ou -- -t AN,AC,AF | \
    bcftools view  --exclude 'AC=0 || F_MISSING>0 || COUNT(GT="het")=N_SAMPLES || COUNT(GT="AA") = N_SAMPLES' | \
    bgzip > $OUT
