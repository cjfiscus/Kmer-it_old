#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=16G
#SBATCH --output=pl%j.stdout
#SBATCH --error=pl%j.stderr
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH --job-name="pl"
#SBATCH -p koeniglab
#SBATCH --array=5-5

# software versions
 samtools 1.8; trimmomatic 0.36; bedtools 2.27.0; jellyfish 2.2.9; bwa 0.7.17

# load required modules (slurm) 
module load trimmomatic/0.36 jellyfish/2.2.9 

##### PARAMETERS ######
# define working directory
WORKINGDIR=./

# path to organellar genomes
ORGANELLAR=/rhome/cfisc004/bigdata/K-mers_Arabidopsis/data/reference/Arabidopsis_thaliana.TAIR10.Mt.Pt.fa

# path to adapter files (for trimmomatic)
ADAPTERS_PE=/rhome/cfisc004/software/Trimmomatic-0.36/adapters/PE_all.fa
ADAPTERS_SE=/rhome/cfisc004/software/Trimmomatic-0.36/adapters/SE_all.fa

# kmer counts will be stored here 
RESULTSDIR=./

# list of sequencing runs
SEQLIST=./sample.txt

##### PIPELINE #####
cd $WORKINGDIR # cd to working directory

# determine sample name
NAME=$(head -n $SLURM_ARRAY_TASK_ID $SEQLIST | tail -n 1 | cut -f1)                    

# define temporary directory
TEMP_DIR=/scratch/cfisc004/$NAME

#-----------------------------------------------------------------------
ftp2=$(head -n $SLURM_ARRAY_TASK_ID $SEQLIST | tail -n 1 | cut -f3) 

LIBTYPE="PE"

if [ "$ftp2" == "" ]; then    							#DEFINE LIBTYPE -> SE OR PE
	LIBTYPE="SE"
fi

FTP1=$(head -n $SLURM_ARRAY_TASK_ID $SEQLIST | tail -n 1 | cut -f2)           # DEFINE FTP1 and FPT2(if exist)  

if [ "$LIBTYPE" == "PE" ]; then                                                 
	FTP2="$ftp2"
fi
#---------------------------------------------------------------------
# make temp directory
mkdir -pv "$TEMP_DIR"
#### ok issue with temp dir here do i needa make a new dir? 

FILE1=$(basename "$FTP1")
FILE2=$(basename "$FTP2") #note to self do i need a if to check if FTP2 exist or nah? crash or nah

if [ "$LIBTYPE" == "PE" ]; then 
	wget -nv "$FTP2"
fi
wget -nv "$FTP1"                                                        # download ftp1 and ftp2 / ALSO locate the filename! via File#

#-----------------------------------------------------------------------
#### TO DO #####
#TESTING WITH TEMP DIR???? HOW 2
#how will i know if the timming work? modify or does it make a new file? 

if [ "$LIBTYPE" == "PE" ]; then

	# Quality/Adapter trimming with Trimmomatic 
	java -jar $TRIMMOMATIC PE -threads 8 \
	$TEMP_DIR/"$FILE1" $TEMP_DIR/"$FILE2" \
	$TEMP_DIR/"$NAME"_1_trimmed_paired.fq.gz $TEMP_DIR/"$NAME"_1_unpaired.fq.gz \
	$TEMP_DIR/"$NAME"_2_trimmed_paired.fq.gz $TEMP_DIR/"$NAME"_2_unpaired.fq.gz \
	ILLUMINACLIP:"$ADAPTERS_PE":2:30:10 \
	LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:24

	# map to mitochondria and plastid genomes (cat'd together) 
	bwa mem -t 8 -M $ORGANELLAR $TEMP_DIR/"$NAME"_1_trimmed_paired.fq.gz $TEMP_DIR/"$NAME"_2_trimmed_paired.fq.gz > $TEMP_DIR/"$NAME".sam

else # single end 
	# Download file
####	wget -nv -P $TEMP_DIR "$FILE" 

	# Quality/Adapter trimming 
	java -jar $TRIMMOMATIC SE -threads 8 \
	$TEMP_DIR/"$FILE1" $TEMP_DIR/"$NAME"_trimmed.fq.gz \
	ILLUMINACLIP:"$ADAPTERS_SE":2:30:10 \
	LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:24
	
	# map to mitochondria and plastid genomes (cat'd together)
	bwa mem -t 8 -M $ORGANELLAR $TEMP_DIR/"$NAME"_trimmed.fq.gz  > $TEMP_DIR/"$NAME".sam

fi

# mapping stats
	#samtools flagstat $TEMP_DIR/"$NAME".sam > $RESULTSDIR/"$NAME"_mapstats.txt

# sam to sorted bam
#samtools view -bS $TEMP_DIR/"$NAME".sam | samtools sort -T $TEMP_DIR/temp_Pt - -o $TEMP_DIR/"$NAME".bam

# extract unmapped reads
#samtools view -f4 -b $TEMP_DIR/"$NAME".bam > $TEMP_DIR/"$NAME".unmapped.bam

# export unmapped reads from original reads 
#bedtools bamtofastq -i $TEMP_DIR/"$NAME".unmapped.bam -fq $TEMP_DIR/"$NAME".unmapped.fq

# Count 12-mers
#jellyfish count -C -m 12 -s 3G -t 8 -o $TEMP_DIR/"$NAME".jf $TEMP_DIR/"$NAME".unmapped.fq
#jellyfish dump -tc $TEMP_DIR/"$NAME".jf > $RESULTSDIR/"$NAME"_12.txt

# clean up
#rm -r $TEMP_DIR
