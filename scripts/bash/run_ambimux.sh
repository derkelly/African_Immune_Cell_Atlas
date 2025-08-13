#!/bin/bash

atac_bam=$1
gex_bam=$2
vcf=$3
gtf=$4
peaks=$5
exclude=$6
samples=$7
bcsd=$8
out=$9

~/tools/ambimux/ambimux \
    --atac-bam "$atac_bam" \
    --rna-bam "$gex_bam" \
    --vcf "$vcf" \
    --gtf "$gtf" \
    --peaks "$peaks" \
    --exclude "$exclude" \
    --samples "$samples" \
    --bc-wl "$bcsd" \
    --rna-mapq 30 \
    --atac-mapq 30 \
    --tx-basic \
    --max-iter 100 \
    --threads 8 \
    --verbose \
    --out "$out" \
    --out-min 100
