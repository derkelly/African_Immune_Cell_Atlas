#!/bin/bash

VIREO_OUTDIR=$1
VCF=$2
N=$3

if [ -e "$VCF" ]; then

    singularity exec --bind /project/tishkofflab ~/tools/Demuxafy.sif vireo \
	-c $VIREO_OUTDIR \
	-d $VCF \
	-o $VIREO_OUTDIR \
	-t GT \
	-N $N \
	--callAmbientRNAs

else
    
    singularity exec --bind /project/tishkofflab ~/tools/Demuxafy.sif vireo \
	-c $VIREO_OUTDIR \
	-o $VIREO_OUTDIR \
	-N $N \
	--callAmbientRNAs

fi
