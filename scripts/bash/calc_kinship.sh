#!/bin/bash

VCF1=$1
VCF1C=$2
VCF2=$3
MERGE=$4
OUT=$5

conda activate muxlet

# test if the quality-filtered file exists
if [ ! -f $VCF1C ]
then
    
    bcftools view  -i  'MIN(FMT/GQ)>10' $VCF1 | bgzip > $VCF1C

    tabix -p vcf $VCF1C
fi

# check if the merged file exists
if [ ! -f $MERGE ]
then
 
    bcftools merge $VCF1C $VCF2 | \
	bcftools filter -e 'GT~"\."' -Oz -o $MERGE
fi

# calculate kinship
plink2 --vcf $MERGE \
    --make-king-table \
    --king-table-filter 0.02 \
    --out $OUT

# get the best guess for each cluster
grep -v ^CLUST $OUT.kin0 | \
    grep CLUST | \
    sort -k2,2 -k6,6nr | \
    awk 'word!=$2{count=1;word=$2} count<=1{print; count++}' > $OUT.kin.best_guess.txt


