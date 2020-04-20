# Kmer-it
A pipeline for counting K-mers from high-throughput sequencing reads.  

This pipeline is designed to count k-mers from reads originating from nuclear genomes. Reads originating from organellar genomes are filtered by mapping all reads to supplied organellar genome(s) and only counting k-mers in unmapped reads.  

Pipeline steps: 
1. Download sequencing reads. 
2. MD5SUM check for file integrity. 
3. Trim reads with trimmomatic. 
4. Map to organellar genome(s), if provided.
5. Extract unmapped reads (if organelle genome is provided)
6. Count K-mers with jellyfish. 

## Dependencies
#### [hstlib/samtools 1.9](https://github.com/samtools/samtools)
#### [trimmomatic 0.36](http://www.usadellab.org/cms/index.php?page=trimmomatic)
#### [bedtools 2.29.2](https://github.com/arq5x/bedtools2)
#### [jellyfish 2.2.9](https://github.com/gmarcais/Jellyfish)
#### [bwa 0.7.17](https://github.com/lh3/bwa)
#### wget (must be installed manually on MacOS, can be done with homebrew or similar)

## Setup 
1. Clone this github repository. 
2. Ensure that dependencies are installed (see above). 
3. Run dl_org_gens.sh to download RefSeq Organellar Genomes from NCBI Organellar Genome Database.
```
sh scripts/dl_org_gens.sh dir/org_gens.fa
```

## Inputs
### Parameter file
Kmer-it requires a parameter file specifying runtime options.  

#### Variables 
SEQ_LIST: Path to list of sequencing runs to process through the pipeline (see Sample file section below).    
THREADS: Number of threads to use for multi-threaded steps (must be integer).  
TEMP_DIR: Path to a temporary directory.  
RM_TEMP_DIR: Set to yes to remove temporary directory at conclusion of pipeline.   
OUT_DIR: Directory where results will be stored.  
RUN_TRIM: Set to yes to run trimmomatic step.  
ADAPTERSPE: Path to file containing paired-end adapters.   
ADAPTERSSE: Path to file containing single-end adapters.  
O_GENOME: Path to organellar genome (mitochondria and/or plastid genomes). Leave blank to skip mapping to organellar genome.  
K: Kmer to count (must be integer).  

### Sample file
Kmer-it requires a file describing the samples to process. Each sequencing run is described on one line. The file must be tab-delimited and contain the following columns (in order):
ID	SAMPLE	FTP	MD5

#### Description of fields
ID: unique id for each sequencing run.   
SAMPLE: sample name. Must be identical for counts from multiple runs to be aggregated.   
FTP: FTP link(s) for file download. Paired-end reads must be contained in separate files and are separated by ";".   
MD5: MD5SUMS for file(s) described in FTP field. If no MD5SUMS are provided this step will be skipped.  

## Usage
1. Count K-mers
```
sh scripts/kmerit.sh params #
```
where params is the parameter file and # is the line to process in the sample file (starts at 2 since 1 is the header).  

2. Aggregate multiple runs (skip if each sample is only sequenced with one run). 
```
sh scripts/combine.sh
```

3. Combine K-mer counts into table. Set DIR to folder where individual K-mer count files are stored. 
```
sbatch 2_mk_count_table.sh ./DIR ./OUT.TBL
```

4. Annotate table with GC content of K-mers in first column.  Run in python3. 
```
RAWCOUNTS=INPUT.TBL
ANNOTCOUNTS=OUTPUT.TBL
paste <(python 3_annot_mers_gc.py <(cut -f1 "$RAWCOUNTS") | cut -f1) "$RAWCOUNTS" > "$ANNOTCOUNTS" 
```

5. Do GC + TMM count normalization. 
```
sbatch 4_normalize.sh ./IN.TBL ./OUT.TBL
```

## Outputs
The following files are created:

## TODO
[ ] Remove temp directory switch
[ ] Combine gc annotation + normalization step
