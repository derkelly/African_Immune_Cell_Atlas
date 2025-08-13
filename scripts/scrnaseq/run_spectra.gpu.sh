#!/bin/bash

IN=$1
OUT=$2
MODEL=$3

module load CUDA/11.7.1

conda activate ctypist

python3 /project/tishkofflab/projects/SingleCellRNA/scripts/scrnaseq/run_spectra.py $IN $OUT $MODEL
