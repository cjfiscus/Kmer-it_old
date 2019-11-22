#!/bin/bash
# Download and format Organellar Genomes from NCBI Organelle Genome Resources
# cjfiscus
# 2019-11-22

# dir/filename to place OrgDB
FILE=$1

# download Mt db
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.1.genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.1.genomic.fna.gz

# download Pt db
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.1.1.genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.2.1.genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.3.1.genomic.fna.gz

# combine files
zcat mitochondrion.1.1.genomic.fna.gz mitochondrion.2.1.genomic.fna.gz \
	plastid.1.1.genomic.fna.gz plastid.2.1.genomic.fna.gz plastid.3.1.genomic.fna.gz \
	> $FILE

# cleanup
rm mitochondrion.*.*.genomic.fna.gz
rm plastid.*.*.genomic.fna.gz

# index with bwa
bwa index $FILE
