#!/bin/bash

IN=$1
OUT=$2

ssh derkelly@consign.pmacs.upenn.edu

cd /archivetape/snutils/
./snretrieve -a $IN

rsync -rplot \
    --inplace \
    --no-partial \
    --whole-file \
    --no-checksum \
    --stats \
    $IN $OUT

exit


