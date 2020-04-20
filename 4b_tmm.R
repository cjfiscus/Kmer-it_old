#!/usr/bin/env Rscript

# install/ load required libs 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, edgeR)

args = commandArgs(trailingOnly=TRUE)

# tmm normalize gc corrected kmer counts

# read in K-mer counts 
counts<-fread(args[1])

# what the data looks like 
dim(counts) 
counts[1:5,1:5]

info<-counts[,1:2]

# TMM normalize counts
factors<-calcNormFactors(DGEList(counts[,3:ncol(counts)]), method="TMM")
rm(counts) # the purge
counts<-cpm(factors)
rm(factors)
out<-as.data.table(cbind(info, counts))
rm(counts)

## write out normalized counts
fwrite(out, file=args[2], sep="\t", quote=F)
