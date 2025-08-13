#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
# default output file
args[2] = "soupx/out"
}

in_d = args[1]
out_d = args[2]

library(Seurat)
library(SoupX)
library(DropletUtils)


sc = load10X(in_d)
sc = autoEstCont(sc)
out = adjustCounts(sc)

# write data
DropletUtils:::write10xCounts(out_d, out)
