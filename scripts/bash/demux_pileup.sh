#!/bin/bash

conda activate muxlet

PATH=~/miniconda3/envs/muxlet/bin:$PATH

IBAM=$1
BCODES=$2
VCF=$3
OBAM=$4
OUT=$5


# keep reads with SNPs and filtered cells from cell ranger
~/tools/popscle_helper_tools/filter_bam_file_for_popscle_dsc_pileup.sh \
    $IBAM \
    $BCODES \
    $VCF \
    $OBAM

# # Running dsc-pileup to generate pileup files; use memory 6* bam
# popscle dsc-pileup \
#     --sam $OBAM \
#     --vcf <(zcat $VCF) \
#     --out $OUT

# dsc-pileup using the Demuxafy singularity image
singularity exec --bind /project/tishkofflab ~/tools/Demuxafy.sif popscle_pileup.py \
    --sam $OBAM \
    --vcf $VCF \
    --group-list $BCODES \
    --out $OUT
