#!/bin/bash -l

# Kmer-it
# cjfiscus
# 2019-03-25

##### PIPELINE #####
# read in arguments from params file
source ./params

# determine sample name
NAME=$(head -n "$1" "$SEQ_LIST" | tail -n 1 | cut -f1)

# determine if seq run is SE or PE
FILE=$(head -n "$1" "$SEQ_LIST" | tail -n 1 | cut -f3)

if [[ "$FILE" == *";"* ]] ; then
        LIBTYPE="PE" # paired end
else
        LIBTYPE="SE" # single end
fi

# work in temp directory
TEMP_DIR="$TEMP_DIR"/"$NAME"
mkdir "$TEMP_DIR"
cd "$TEMP_DIR"

if [ $LIBTYPE == "PE" ]
then # paired end 
        # Download files
	INDEX=1 # 1 is forward, 2 is reverse
	for i in $(echo $FILE | tr ";" "\n")
	do	
		echo "downloading" "$i" 
		wget "$i" -O "$NAME"_"$INDEX".fastq.gz
		INDEX=$((INDEX + 1))
	done

	# check MD5sums
	echo "verifying checksums..."
	INDEX=1
	MD5SUMS=$(head -n "$1" "$SEQLIST" | tail -n 1 | cut -f4)
	if [ -z "$MD5SUMS" ]
	then 
		echo "skipping MD5sum check..."
	else

		for i in $(echo "$MD5SUMS" | tr ";" "\n")
		do 
			echo "$i" "$NAME"_"$INDEX".fastq.gz >> chk.md5
    			INDEX=$((INDEX + 1))
		done	

		if md5sum --status -c chk.md5; then
			# continue		
			echo "SUMS GOOD"
		else
			# stop script
			echo "SUMS BAD"
			exit 1
		fi
	fi

	if [[ $RUN_TRIM = "yes" ]]
	then 
		# Quality/Adapter trimming
		echo "trimming with trimmomatic..."
		java -jar $TRIMMOMATIC PE -threads "$THREADS" \
		"$NAME"_1.fastq.gz "$NAME"_2.fastq.gz \
		"$NAME"_1_trimmed_paired.fq.gz "$NAME"_1_unpaired.fq.gz \
		"$NAME"_2_trimmed_paired.fq.gz "$NAME"_2_unpaired.fq.gz \
		ILLUMINACLIP:"$ADAPTERSPE":2:30:10 \
		LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:36

	else
		mv "$NAME"_1.fastq.gz "$NAME"_1_trimmed_paired.fq.gz
		mv "$NAME"_2.fastq.gz "$NAME"_2_trimmed_paired.fq.gz
	fi

	if [ -z "$REF_GENOME" ]
	then 
		echo "No mapping to reference genome"

	else
		# map to reference genome
		echo "mapping with bwa..."
		bwa mem -t 8 -M $REF_GENOME "$NAME"_1_trimmed_paired.fq.gz \
			"$NAME"_2_trimmed_paired.fq.gz > "$NAME"_gen.sam

	fi

	if [ -z "$O_GENOME" ]
	then 
		echo "No mapping to organellar genome"
	
	else
		# map to organellar genome 
		bwa mem -t 8 -M $O_GENOME "$NAME"_1_trimmed_paired.fq.gz \
			"$NAME"_2_trimmed_paired.fq.gz > $TEMP_DIR/"$NAME"_org.sam
	fi

else # single end 
	# Download file
	echo "downloading" "$FILE"	
	wget "$FILE" -O "$NAME".fastq.gz
	
	# check MD5sum
	MD5SUMS=$(head -n "$1" $SEQLIST | tail -n 1 | cut -f4)
	if [ -z "$MD5SUMS" ]
	then 
		echo "skipping MD5sum check..."
	else
		echo "verifying checksums..."

		echo "$MD5SUMS" "$NAME".fastq.gz >> chk.md5
	
		if md5sum --status -c chk.md5; then
        		# continue
        		echo "SUMS GOOD"
        	else 
        		# stop script
        		echo "SUMS BAD"
        		exit 1
        	fi
	fi 

	if [[ $RUN_TRIM = "yes" ]]
	then 
		# Quality/Adapter trimming
		echo "trimming with trimmomatic..."
		java -jar $TRIMMOMATIC SE -threads 8 \
		"$NAME".fastq.gz "$NAME"_trimmed.fq.gz \
		ILLUMINACLIP:"$ADAPTERSSE":2:30:10 \
		LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:36

	else 
		mv "$NAME".fastq.gz "$NAME"_trimmed.fq.gz
	fi

	if [ -z "$REF_GENOME" ]
	then 
		# map to reference genome
		echo "mapping with bwa..."
		bwa mem -t 8 -M $REF_GENOME "$NAME"_trimmed.fq.gz  > "$NAME"_gen.sam
	fi 

	if [ -z "$O_GENOME" ]
	then 
		# map to organellar genome
		bwa mem -t 8 -M $O_GENOME "$NAME"_trimmed.fq.gz  > "$NAME"_org.sam

	fi 
fi

# sam to sorted bam
#echo "samtools sam to bam"
#samtools view -bS "$NAME"_gen.sam | samtools sort -T temp_Pt - -o "$NAME"_gen.bam
#samtools view -bS "$NAME"_org.sam | samtools sort -T temp_Pt - -o "$NAME"_org.bam

# mapping stats
#echo "samtools flagstat"
#samtools flagstat "$NAME"_gen.bam > $RESULTSDIR/mappingstats/genome/"$NAME"_gen_mapstats.txt
#samtools flagstat "$NAME"_org.bam > $RESULTSDIR/mappingstats/organellar/"$NAME"_org_mapstats.txt

# index bam
#echo "samtools indexing bam"
#samtools index "$NAME"_gen.bam

# calculate coverage of ref per base and in 1kb windows
#echo "calculating coverage with mosdepth..."
#mosdepth -t "$THREADS" -b 1000 "$RESULTSDIR"/coverage/genome/"$NAME" "$NAME"_gen.bam

# subset reads that did not map to organelles
#echo "extracting unmapped reads..."
#samtools view -f4 -b "$NAME"_org.bam > "$NAME".unmapped.bam

# export these unmapped reads
#bedtools bamtofastq -i "$NAME".unmapped.bam -fq "$NAME".unmapped.fq

# Count 10-mers in reads that did not map to organelles 
#echo "counting K-mers with jellyfish"
#jellyfish count -C -m "$K" -s 3G -t "$THREADS" -o "$NAME".jf "$NAME".unmapped.fq 
#jellyfish dump -tc "$NAME".jf > $OUT_DIR/"$NAME".txt
