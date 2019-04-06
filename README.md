# Kmer-it
A pipeline for counting K-mers from high-throughput sequencing reads.  
This pipeline is designed to count k-mers from reads originating from nuclear genomes. Reads originating from organellar genomes are filtered by mapping all reads to supplied organellar genome(s) and only counting k-mers in unmapped reads.  

Pipeline steps: 
1. Download sequencing reads with axel. 
2. MD5SUM check for file integrity. 
3. Trim reads with trimmomatic. 
4. Map to reference genome. 
5. Map to organellar genome(s).
6. Calculate coverage with mosdepth (if reference genome is provided).
7. Extract unmapped reads (if organelle genome is provided)
8. Count K-mers with jellyfish. 
9. Repeat assembly with REPdenovo (reference genome must be provided). 

## Dependencies
#### [axel 2.16.1](https://github.com/axel-download-accelerator/axel)
#### [hstlib/samtools 1.9](https://github.com/samtools/samtools)
#### [trimmomatic 0.36](http://www.usadellab.org/cms/index.php?page=trimmomatic)
#### [bedtools 2.28.0](https://github.com/arq5x/bedtools2)
#### [jellyfish 2.2.29](https://github.com/gmarcais/Jellyfish)
#### [bwa 0.7.17](https://github.com/lh3/bwa)
#### [REPdenovo](https://github.com/Reedwarbler/REPdenovo)

## Inputs
### Parameter file
Kmer-it requires a parameter file specifying runtime options.  

#### Variables 
SEQ_LIST: Path to list of sequencing runs to process through the pipeline (see Sample file section below).  
THREADS: Number of threads to use for multi-threaded steps (must be integer). 
TEMP_DIR: Path to a temporary directory
RM_TEMP_DIR: Set to yes to remove temporary directory at conclusion of pipeline. 
OUT_DIR: Directory where results will be stored. 
RUN_TRIM: Set to yes to run trimmomatic step. 
ADAPTERSPE: Path to file containing paired-end adapters. 
ADAPTERSSE: Path to file containing single-end adapters. 
REF_GENOME: Path to reference genome. Leave blank to skip mapping to reference genome. 
O_GENOME: Path to organellar genome (mitochondria and/or plastid genomes). Leave blank to skip mapping to organellar genome. 
K: Kmer to count (must be integer). 
REP_ASSEM: Set to yes to do repeat assembly with REPdenovo. 
REPDENOVO: Path to REPdenovo script (main.py). 
REP_CONFIG: Path to REPdenovo configuration file. 

### Sample file
Kmer-it requires a file describing the samples to process. Each sequencing run is described on one line. The file must be tab-delimited and contain the following columns (in order):
ID	SAMPLE	FTP	MD5

#### Description of fields
ID: unique id for each sequencing run.   
SAMPLE: sample name. Must be identical for counts from multiple runs to be aggregated.   
FTP: FTP link(s) for file download. Paired-end reads must be contained in separate files and are separated by ";".   
MD5: MD5SUMS for file(s) described in FTP field. If no MD5SUMS are provided this step will be skipped.  

## Usage
Running Kmer-it consists of two steps. 

1. Count K-mers
```
sh kmerit.sh params #
```
where params is the parameter file and # is the line to process in the sample file (starts at 2 since 1 is the header).  

2. Aggregate multiple runs (skip if each sample is only sequenced with one run). 
```
sh combine.sh
```
## Outputs
The following files are created:

## TODO
[ ] Remove temp directory switch
[ ] Test repeat assembly
