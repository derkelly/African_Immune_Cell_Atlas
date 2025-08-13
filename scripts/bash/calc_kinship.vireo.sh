#!/bin/bash

VCF1=$1
VCF2=$2
MERGE=$3
OUT=$4

conda activate muxlet

# check if the merged file exists
if [ ! -f $MERGE ]
then

    bcftools merge --force-samples $VCF1 $VCF2 | \
    	bcftools filter -e 'GT~"\."' | bgzip > $MERGE

else

    FSIZE=$(stat -c%s "$MERGE")

    # if the file is small then re-run merge
    if [[ $FSIZE < 1000000 ]]
    then

	bcftools merge --force-samples $VCF1 $VCF2 | \
    	    bcftools filter -e 'GT~"\."' | bgzip > $MERGE

    fi
fi

# calculate kinship
plink2 --vcf $MERGE \
    --make-king-table \
    --king-table-filter 0.02 \
    --out $OUT

# get the best guess for each cluster
grep -v ^donor $OUT.kin0 | \
    grep "donor" | \
    sort -k2,2 -k6,6nr | \
    awk 'word!=$2{count=1;word=$2} count<=1{print; count++}' > $OUT.kin.best_guess.txt
