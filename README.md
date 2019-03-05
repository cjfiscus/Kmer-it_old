# Kmer-it
A pipeline for counting K-mers from high-throughput sequencing reads. Currently configured to run on CentOS and slurm queueing system. 

## Dependencies
#### wget
#### samtools 1.8
#### trimmomatic 0.36
#### bedtools 2.27.0
#### jellyfish 2.2.29
#### bwa 0.7.17

## Usage

1_pipeline.sh samples.txt

## TODO
 [ ] Add additional column in test file for identifier  
 [ ] Check that all identifers are unique  
 [ ] Add in NA check for FTP2  
 [ ] Add assembly step  
 [ ] Report % aligning to genome  
 [ ] Move all options to a params file  
 [ ] Write second script to combine samples  
