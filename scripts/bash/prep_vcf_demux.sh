#!/bin/bash

IDS=$1
GENOS=$2
OUT=$3

bcftools view \
    --samples-file $IDS \
    -Ou \
    $GENOS | \
    bcftools +fill-tags -Ou -- -t AN,AC,AF | \
    bcftools annotate -x FILTER,INFO/MAF,INFO/ExcHet,INFO/F_MISSING -Ou | \
    bcftools view --exclude 'AC=0 || F_MISSING>0 || COUNT(GT="het")=N_SAMPLES || COUNT(GT="AA") = N_SAMPLES' | \
    bgzip > $OUT
