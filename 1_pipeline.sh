#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=40G
#SBATCH --output=std/pl%j.stdout
#SBATCH --error=std/pl%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="pl"
#SBATCH -p batch
#SBATCH --array=1-15

# software versions
# samtools 1.8; trimmomatic 0.36; bedtools 2.27.0; jellyfish 2.2.9; bwa 0.7.17

# load required modules (slurm) 
module load trimmomatic/0.36 jellyfish/2.2.9 

##### PARAMETERS ######
# define working directory
WORKINGDIR=../../

# path to organellar genomes
ORGANELLAR=/rhome/cfisc004/bigdata/K-mers_Arabidopsis/data/reference/Arabidopsis_thaliana.TAIR10.Mt.Pt.fa

# path to adapter files (for trimmomatic)
ADAPTERS_PE=/rhome/cfisc004/software/Trimmomatic-0.36/adapters/PE_all.fa
ADAPTERS_SE=/rhome/cfisc004/software/Trimmomatic-0.36/adapters/SE_all.fa

# kmer counts will be stored here 
RESULTSDIR=/rhome/cfisc004/bigdata/ma_lines/results/

# list of sequencing runs
SEQLIST=/rhome/cfisc004/bigdata/ma_lines/data/MA_lines.txt

##### #####

##### PIPELINE #####
cd $WORKINGDIR # cd to working directory

# get filenames from list 
FILE=$(head -n $SLURM_ARRAY_TASK_ID $SEQLIST | tail -n 1 | cut -f2)

# determine filenames
# FILE1 is forward read if paired end, FILE2 is reverse read
for i in $(seq 1 $(echo $FILE | tr ";" "\n" | wc -l))
do
        declare "FILE$i"=$(basename $(echo $FILE | tr ";" "\n" | head -n "$i" | tail -n 1))
done

# determine sample name
NAME=$(head -n $SLURM_ARRAY_TASK_ID $SEQLIST | tail -n 1 | cut -f1)

# define temporary directory
TEMP_DIR=/scratch/cfisc004/$NAME

# check if library is single end or paired end  
if [[ "$FILE" == *";"* ]] ; then
	LIBTYPE="PE" # paired end 
else
	LIBTYPE="SE" # single end 
fi

# make temp directory
mkdir -pv "$TEMP_DIR"

if [ $LIBTYPE == "PE" ]
then # paired end 
	# Download files
        for i in $(echo $FILE | tr ";" "\n")
	do
		wget -nv -P $TEMP_DIR "$i"
	done

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
	wget -nv -P $TEMP_DIR "$FILE" 

	# Quality/Adapter trimming 
	java -jar $TRIMMOMATIC SE -threads 8 \
	$TEMP_DIR/"$FILE1" $TEMP_DIR/"$NAME"_trimmed.fq.gz \
	ILLUMINACLIP:"$ADAPTERS_SE":2:30:10 \
	LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:24
	
	# map to mitochondria and plastid genomes (cat'd together)
	bwa mem -t 8 -M $ORGANELLAR $TEMP_DIR/"$NAME"_trimmed.fq.gz  > $TEMP_DIR/"$NAME".sam

fi

# mapping stats
samtools flagstat $TEMP_DIR/"$NAME".sam > $RESULTSDIR/"$NAME"_mapstats.txt

# sam to sorted bam
samtools view -bS $TEMP_DIR/"$NAME".sam | samtools sort -T $TEMP_DIR/temp_Pt - -o $TEMP_DIR/"$NAME".bam

# extract unmapped reads
samtools view -f4 -b $TEMP_DIR/"$NAME".bam > $TEMP_DIR/"$NAME".unmapped.bam

# export unmapped reads from original reads 
bedtools bamtofastq -i $TEMP_DIR/"$NAME".unmapped.bam -fq $TEMP_DIR/"$NAME".unmapped.fq

# Count 12-mers
jellyfish count -C -m 12 -s 3G -t 8 -o $TEMP_DIR/"$NAME".jf $TEMP_DIR/"$NAME".unmapped.fq
jellyfish dump -tc $TEMP_DIR/"$NAME".jf > $RESULTSDIR/"$NAME"_12.txt

# clean up
rm -r $TEMP_DIR
