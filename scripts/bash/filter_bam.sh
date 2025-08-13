#!/bin/bash

BAM=$1
BCODES=$2
VCF=$3
OUT=$4


# keep reads with SNPs and filtered cells from cell ranger
~/tools/popscle_helper_tools/filter_bam_file_for_popscle_dsc_pileup.sh \
    $BAM \
    $BCODES \
    $VCF \
    $OUT

