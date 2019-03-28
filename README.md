# Kmer-it
A pipeline for counting K-mers from high-throughput sequencing reads. 

## Dependencies
#### [hstlib/samtools 1.9](https://github.com/samtools/samtools)
#### [trimmomatic 0.36](http://www.usadellab.org/cms/index.php?page=trimmomatic)
#### [bedtools 2.28.0](https://github.com/arq5x/bedtools2)
#### [jellyfish 2.2.29](https://github.com/gmarcais/Jellyfish)
#### [bwa 0.7.17](https://github.com/lh3/bwa)

## Usage

1_pipeline.sh samples.txt

## TODO
 [ ] Add additional column in test file for identifier  
 [ ] Add in NA check for FTP2  
 [ ] Add assembly step  
 [ ] Report % aligning to genome  
 [ ] Write second script to combine samples  
