#!/bin/bash

IN=$1
OUT=$2

ssh derkelly@mercury.pmacs.upenn.edu

rsync -rplot \
    --inplace \
    --no-partial \
    --whole-file \
    --no-checksum \
    --stats \
    $IN $OUT

exit
